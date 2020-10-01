<html>
<?php 
if(file_exists("/update.log")) { 
	$zawartosc_logu = file_get_contents('/update.log');
	$czy_dziala = shell_exec('ps -ef | grep -v grep | grep OmicSelector-updater | wc -l');
	$skonczone = 0; if (strpos($zawartosc_logu, 'update is finished') !== false) { $skonczone = 1; } } else { $msg = urlencode("The update was not initialized. Please run it again."); header("Location: /?msg=" . $msg); die(); }
?>

<head>
    <title>OmicSelector (updater)</title>
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
    $( document ).ready(function() {
		$('html, body').scrollTop($(document).height());
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
    <p>Update log:
    <pre><?php echo $zawartosc_logu; ?></pre></p>
    
    <?php if($skonczone == 1) { ?>
	<p id="msg"><b>The update is finished.</b> Please go back to the app. If you have any active notebooks running, you may need to restart kernel for new features.<br />More details on new features: <a href="https://github.com/kstawiski/OmicSelector" target="_blank">https://github.com/kstawiski/OmicSelector</a>.</p>
    <p></p><a href="/" onclick="waitingDialog.show('Going back...');" class="btn btn-success"><i class="fas fa-undo"></i>&emsp;Go back</a></p>
	<?php } else { ?>
		<p id="msg"><b>The update is still in progress...</b> Please do not use the app and don't leave this page. Running status: <code><?php echo $czy_dziala; ?></code></p>
		<meta http-equiv="refresh" content="3">
	<?php } ?>
    <p>&emsp;</p>
	</div>
  </body>
</html>