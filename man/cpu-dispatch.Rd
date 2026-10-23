\name{cpu-dispatch}
\alias{cpu-dispatch}
\title{Foundation-model CPU dispatch and provenance}
\description{OmicSelector sets portable x86-64 CPU dispatch defaults when the
package loads, before importing Python backends. Explicit environment settings
are preserved, including empty values.}
\details{
The defaults are \code{ATEN_CPU_CAPABILITY=avx2},
\code{MKL_ENABLE_INSTRUCTIONS=AVX2}, \code{ONEDNN_MAX_CPU_ISA=AVX2},
and \code{DNNL_MAX_CPU_ISA=AVX2}. On Linux hosts whose \file{/proc/cpuinfo}
reports AVX-512, \code{OPENBLAS_CORETYPE=Zen} selects the validated AVX2
OpenBLAS implementation. Other hosts retain native OpenBLAS dispatch.
They require an AVX2-capable x86-64 CPU. Other architectures retain their native
dispatch. Library versions, thread counts and other numerical settings still
need to match for reproducible scores.

Models fitted by \code{fit_tabpfn}, \code{fit_tabdpt} and \code{fit_tabicl}
store \code{torch.backends.cpu.get_cpu_capability()} and these environment
settings in \code{model$provenance$cpu_dispatch}. Scoring compares the active
capability and settings with that record and warns on a mismatch.
Use \code{options(OmicSelector.cpu_dispatch_strict = TRUE)} to reject mismatches.
Strict mode also rejects legacy models without dispatch provenance; the default
allows these models because their fitted capability is unknown.

If another package initialized torch or a numerical library first, changing
environment variables may be too late. Restart R and set dispatch before the
first import. The provenance records the active torch capability; it does not
claim to introspect the initialized ISA of every BLAS or oneDNN kernel.
No probability correction, rounding or score calibration is performed.
}
\seealso{\code{\link{fit_tabpfn}}, \code{\link{fit_tabdpt}},
\code{\link{fit_tabicl}}}
