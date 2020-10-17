<html>

<head>
    <title>OmicSelector</title>
    <script src="https://code.jquery.com/jquery-3.5.1.js"
   
  crossorigin="anonymous"></script>
    <script src="jquery-ui.js" type="text/javascript"></script>
    <link rel="stylesheet" href="bootstrap.min.css"
          />
    <link rel="stylesheet" href="bootstrap-theme.min.css"
          crossorigin="anonymous" />
    <script src="bootstrap.min.js"
         >
    </script>
    <meta charset="utf-8" />
    <meta http-equiv="X-UA-Compatible" content="IE=edge" />
    <meta name="viewport" content="width=device-width, initial-scale=1" />
    <meta name="description" content="OmicSelector - a tool for selecting great miRNA biomarkers." />
    <meta name="author" content="Konrad Stawiski (konrad@konsta.com.pl)" />
    <link rel="stylesheet" href="css/starter-template.css" />
    <script src="https://cdnjs.cloudflare.com/ajax/libs/font-awesome/5.14.0/js/all.min.js"
         ></script>
    <!-- Global site tag (gtag.js) - Google Analytics -->
    <script async src="https://www.googletagmanager.com/gtag/js?id=UA-53584749-8"></script>
    <script>
    window.dataLayer = window.dataLayer || [];
    function gtag(){dataLayer.push(arguments);}
    gtag('js', new Date());

    gtag('config', 'UA-53584749-8');
    </script>
    <script type="text/javascript">
    $(".btn-success").click(function (event) {
        waitingDialog.show('Processing.. Please wait...');
            });

    var waitingDialog = waitingDialog || (function ($) { 'use strict';

	// Creating modal dialog's DOM
	var $dialog = $(
		'<div class="modal fade" data-backdrop="static" data-keyboard="false" tabindex="-1" role="dialog" aria-hidden="true" style="padding-top:15%; overflow-y:visible;">' +
		'<div class="modal-dialog modal-m">' +
		'<div class="modal-content">' +
			'<div class="modal-header"><h3 style="margin:0;"></h3></div>' +
			'<div class="modal-body">' +
				'<div class="progress progress-striped active" style="margin-bottom:0;"><div class="progress-bar" style="width: 100%"></div></div>' +
			'</div>' +
		'</div></div></div>');

	return {
		/**
		 * Opens our dialog
		 * @param message Custom message
		 * @param options Custom options:
		 * 				  options.dialogSize - bootstrap postfix for dialog size, e.g. "sm", "m";
		 * 				  options.progressType - bootstrap postfix for progress bar type, e.g. "success", "warning".
		 */
		show: function (message, options) {
			// Assigning defaults
			if (typeof options === 'undefined') {
				options = {};
			}
			if (typeof message === 'undefined') {
				message = 'Loading';
			}
			var settings = $.extend({
				dialogSize: 'm',
				progressType: '',
				onHide: null // This callback runs after the dialog was hidden
			}, options);

			// Configuring dialog
			$dialog.find('.modal-dialog').attr('class', 'modal-dialog').addClass('modal-' + settings.dialogSize);
			$dialog.find('.progress-bar').attr('class', 'progress-bar');
			if (settings.progressType) {
				$dialog.find('.progress-bar').addClass('progress-bar-' + settings.progressType);
			}
			$dialog.find('h3').text(message);
			// Adding callbacks
			if (typeof settings.onHide === 'function') {
				$dialog.off('hidden.bs.modal').on('hidden.bs.modal', function (e) {
					settings.onHide.call($dialog);
				});
			}
			// Opening dialog
			$dialog.modal();
		},
		/**
		 * Closes dialog
		 */
		hide: function () {
			$dialog.modal('hide');
		}
	};

})(jQuery);
    </script>
</head>

<body>
    <div class="container">
        <div class="starter-template">
            <p>
                <center><img src="logo.png" width="70%" />
            </p>
            <p><br></p>
        </div>
        <p>Welcome to <b>OmicSelector</b>. This form can be used to set up the new analysis.</p>
        <div class="panel-group">
            <?php if ($_GET["msg"] != "") { ?>
            <div class="panel panel-danger">
                <div class="panel-heading"><i class="fas fa-exclamation-triangle"></i></i>&emsp;&emsp;MESSAGE</div>
                <div class="panel-body"><b><?php echo htmlentities($_GET['msg']); ?></b></div>
            </div>
            <?php } ?>

            <div class="panel panel-success">
                <div class="panel-heading"><i class="fas fa-cloud-upload-alt"></i>&emsp;&emsp;Start new analysis!</div>
                <div class="panel-body">

                    <p><b>How to prepare files for the analysis?</b>
                    <ul>
                        <li>All variables which (column) name starts with <code>hsa</code> will be considered features of interest. Features of interest must be numeric and have no missing values.</li>
                        <li>All variables which name doesn't start with <code>hsa</code> will be considered as metadata.</li>
                        <li>The dataset has to have the variable named <code>Class</code> with no missing data and with values expicitly encoded as <code>Cancer</code> or <code>Control</code> cases. If your problem does not involve cancer detecton, please temporary change your encoding and interpret the results accordingly.</li>
                        <li>The pipeline will perform random data spltting in ratio 60% (training set) : 20% (testing set) : 20% (validation set). If you wish to enforce your way of spltting, the submitted file should have the variable named <code>mix</code> with values <code>train</code>, <code>test</code> and <code>valid</code>.</li>
                    </ul>
                    </p>
                    <p><a href="https://github.com/kstawiski/OmicSelector/blob/master/example/Elias2017.csv" target="_blank">See exemplary file: <code>Elias2017.csv</code></a>. This file originates from our paper <a href="https://elifesciences.org/articles/28932" target="_blank">Elias et al. 2017</a>.
                    </p>
                    <hr>
                    <p><b>Upload the file for the analysis:</b></p>
                    <form action="process.php?type=new_analysis" method="post" enctype="multipart/form-data">
                        <p>Type of expression data:
                        <select name="type" id="type">
                            <option value="logtpm">Log10-transformed values (e.g. logTPM from NGS experiments)</option>
                            <option value="deltact">Crude values (e.g. deltaCt values from qPCR analysis)</option>
                        </select>
                        </p>

                        <p>Select <code>.csv</code> (comma-seperated values) or <code>.xlsx</code> (Microsoft Excel file) file to start the analysis:</p>
                        <input type="file" class="form-control-file" id="fileToUpload" name="fileToUpload"><br />

                        <button type="submit" class="btn btn-success" value="Upload" name="submit" onclick="waitingDialog.show('Uploading and performing initial check...');">
                        <i class="fas fa-upload"></i>&emsp;Upload
                        </button>&emsp;<a href="/" onclick="waitingDialog.show('Going back...');" class="btn btn-success"><i class="fas fa-sign-out-alt"></i>&emsp;Exit</a>
                    </form>

                </div>
            </div>

    </div>

    <hr />
    <footer class="footer">
        <div class="container">
            <span class="text-muted">OmicSelector by Konrad Stawiski and Marcin Kaszkowiak&emsp;&emsp;&emsp;&emsp;<i
                    class="fas fa-envelope"></i> konrad@konsta.com.pl&emsp;&emsp;&emsp;<i
                    class="fas fa-globe-europe"></i>
                <a href="https://biostat.umed.pl" taret="_blank">https://biostat.umed.pl</a>&emsp;&emsp;&emsp;<i
                    class="fab fa-github"></i> <a href="https://kstawiski.github.io/OmicSelector/"
                    target="_blank">https://kstawiski.github.io/OmicSelector/</a></span>
                    <p>&emsp;</p>
        </div>
    </footer>
    <!-- /.container -->
</body>

</html>


