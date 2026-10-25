"""TabDPT adapter with CPU-fitted preprocessing and a frozen CPU projection."""
import numpy as np
import torch
from tabdpt import TabDPTClassifier
from foundation_inference import cpu_projection


class FrozenTabDPTClassifier(TabDPTClassifier):
    def fit(self, X, y, seed=42):
        # Fit the backend's imputer/scaler on CPU, suppressing its device-local
        # randomized PCA. Only fitted training rows enter the projection.
        requested_device = self.device
        self.device = "cpu"
        self.feature_reduction = "subsample"
        if self.scaler is not None and hasattr(self.scaler, "random_state"):
            self.scaler.set_params(random_state=int(seed))
        super().fit(X, y)
        self.feature_reduction = "pca"
        if self.n_features > self.max_features:
            self.V, diagnostics = cpu_projection(
                self.X_train, self.max_features, seed, return_metadata=True)
            for name, value in diagnostics.items():
                setattr(self, "projection_" + name, value)
        else:
            self.V = None
            for name in ("rank", "q", "sign_min_margin", "rank_tolerance",
                         "min_retained_ratio"):
                if hasattr(self, "projection_" + name):
                    delattr(self, "projection_" + name)
        self.projection_seed = int(seed)
        self.projection_device = "cpu"
        self.inference_precision = "float32"
        self._context_cpu = self._project_cpu(self.X_train)
        self.device = requested_device
        return self

    def _project_cpu(self, X):
        x = torch.as_tensor(np.array(X, copy=True), device="cpu", dtype=torch.float32)
        return x if self.V is None else x @ self.V.cpu()

    def _prepare_prediction(self, X, class_perm=None, seed=None):
        from sklearn.utils.validation import check_is_fitted
        check_is_fitted(self)
        if self.missing_indicators:
            inds = np.isnan(X)[:, self.has_missing_indicator].astype(float)
            X = np.hstack((X, inds))
        self.X_test = self.imputer.transform(X)
        if self.scaler is not None:
            self.X_test = self.scaler.transform(self.X_test)
            if self.normalizer == "quantile-uniform":
                self.X_test = 2 * self.X_test - 1
        train_x = self._context_cpu.to(self.device)
        test_x = self._project_cpu(self.X_test).to(self.device)
        train_y = torch.as_tensor(self.y_train, device=self.device, dtype=torch.float32)
        if class_perm is not None:
            inv_perm = torch.as_tensor(np.argsort(class_perm), device=self.device)
            train_y = inv_perm[train_y.long()].float()
        return train_x, train_y, test_x

    def to(self, device):
        # CPU preprocessing and projection stay canonical across device moves.
        self.device = str(device)
        self.model.to(device)
        return self
