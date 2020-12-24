<?php
exec("chown -R OmicSelector /OmicSelector");
require_once 'class.formr.php';
if (!file_exists('/OmicSelector/var_status.txt')) { file_put_contents('/OmicSelector/var_status.txt', "[0] INITIAL (UNCONFIGURED)"); } // Wyjściowy status.
$status = file_get_contents('/OmicSelector/var_status.txt');
$version = file_get_contents('/version.txt');
$pid = shell_exec("ps -ef | grep -v grep | grep OmicSelector-task | awk '{print $2}'");
if($pid != "") { header("Location: /inprogress.php"); }

session_start();
$prev_analysis = "";
if(isset($_SESSION["analysis_id"]))
{
    $prev_analysis = $_SESSION["analysis_id"];
}


?>


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

    $( "body" ).prepend( '<div id="preloader"><div class="spinner-sm spinner-sm-1" id="status"></div></div>' );
$(window).on('load', function() { // makes sure the whole site is loaded
  $('#status').fadeOut(); // will first fade out the loading animation
  $('#preloader').delay(350).fadeOut('slow'); // will fade out the white DIV that covers the website.
  $('body').delay(350).css({'overflow':'visible'});
})

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

        <div class="panel-group">
            <?php if ($_GET["msg"] != "") { ?>
            <div class="panel panel-danger">
                <div class="panel-heading"><i class="fas fa-exclamation-triangle"></i></i>&emsp;&emsp;MESSAGE</div>
                <div class="panel-body"><b><?php echo htmlentities($_GET['msg']); ?></b></div>
            </div>
            <?php } ?>
            <div class="panel panel-default">
                <div class="panel-heading"><i class="fas fa-info"></i>&emsp;&emsp;Welcome to OmicSelector</div>
                <div class="panel-body"><p>Welcome to <b>OmicSelector</b> - the software intended to find the best biomarker signiture based on NGS or other omic modalities (miRNA-seq, RNA-seq).</p>

        <p>Your current version of software: <code>OmicSelector v1.0.<?php echo $version; ?></code></p>
        <?php
            if (file('/PUBLIC', FILE_IGNORE_NEW_LINES)[0] == "1") {
                echo '<p><b>This is a public (demo) version of software. Some hacking options are restricted.</b></p>';
            }
        ?>
        </div>
            </div>
            <div class="panel panel-primary">
                <div class="panel-heading"><i class="fas fa-chart-pie"></i>&emsp;&emsp;Analysis</div>
                <div class="panel-body">
                <p>You can start new analysis or resume your previous one.</p>
                <table class="table">
                    <tr><td>
                    <h4>Start new analysis:</h4>
                    <p><a href="/start.php" role="button" class="btn btn-primary"><i class="fas fa-plus"></i>&emsp;New OmicSelector analysis</a></p>
                    <p><a href="/start_custom.php" role="button" class="btn btn-warning"><i class="fas fa-plus"></i>&emsp;New custom analysis</a></p>
                    </td>
                    <td>
                        <h4>Resume analysis:</h4>
                <p><form action="/analysis.php" method="get">
                    <p><input type="text" id="id" name="id" placeholder="Provide analysis ID" class="form-control" value="<?php echo $prev_analysis; ?>" autocomplete="analysisid name"></p>
                <p><button type="submit" class="btn btn-success" value="Upload" name="submit" onclick="waitingDialog.show('Loading...');"><i class="fas fa-folder-open"></i>&emsp;Resume analysis</button></p>
                </form></p>
            </td></tr>
            </table>
                </div>
            </div>


            <div class="panel panel-warning">
                <div class="panel-heading"><i class="fas fa-puzzle-piece"></i>&emsp;&emsp;Preprocessing extensions</div>
                <div class="panel-body">Coming soon...</div>
            </div>



        <div class="panel panel-default">
                <div class="panel-heading"><i class="fas fa-bars"></i>&emsp;&emsp;Additional tools</div>
                <div class="panel-body"><button type="button" class="btn btn-info" data-toggle="modal"
                        data-target="#modalYT"><i class="fas fa-tv"></i>&emsp;System monitor</button>&emsp;<button type="button" class="btn btn-info" data-toggle="modal" data-target="#modalYT2"><i class="fas fa-terminal"></i>&emsp;Shell</button>&emsp;
                        <a href="monitor/" target="_blank" role="button" class="btn btn-info"><i class="fas fa-server"></i>&emsp;Hardware</a>&emsp;<a href="process.php?type=init_update" role="button" onclick="waitingDialog.show('Starting update...');" class="btn btn-primary"><i class="fas fa-arrow-up"></i></i>&emsp;Update</a>
                    &emsp;<a href="e/notebooks/OmicSelector/vignettes/Tutorial.Rmd" role="button" class="btn btn-primary" target="_blank"><i class="fas fa-graduation-cap"></i>&emsp;Learn R package</a>
                    <?php
                    // Advanced features in private version:
                    if (file('/PUBLIC', FILE_IGNORE_NEW_LINES)[0] != "1") {
                        echo '&emsp;<a href="e/tree" role="button" class="btn btn-primary" target="_blank"><i class="fas fa-puzzle-piece"></i>&emsp;Advanced features</a>';
                    }
                    ?>
                    </div>
            </div>

    </div>
    <!--Modal: Name-->
    <div class="modal fade" id="modalYT" tabindex="-1" role="dialog" aria-labelledby="myModalLabel" aria-hidden="true">
        <div class="modal-dialog modal-lg" role="document">

            <!--Content-->
            <div class="modal-content">

                <!--Body-->
                <div class="modal-body mb-0 p-0">

                    <div class="embed-responsive embed-responsive-16by9 z-depth-1-half">
                        <iframe class="embed-responsive-item" src="top.php" allowfullscreen></iframe>
                    </div>

                </div>

                <!--Footer-->
                <div class="modal-footer justify-content-center">
                    <span class="mr-4">Running <code>top</code> every 2 seconds...</span>

                    <button type="button" class="btn btn-outline-primary btn-rounded btn-md ml-4"
                        data-dismiss="modal">Close</button>

                </div>

            </div>
            <!--/.Content-->

        </div>
    </div>
    <!--Modal: Name-->


        <!--Modal: Name-->
        <div class="modal fade" id="modalYT2" tabindex="-1" role="dialog" aria-labelledby="myModalLabel" aria-hidden="true">
        <div class="modal-dialog modal-lg" role="document">

            <!--Content-->
            <div class="modal-content">

                <!--Body-->
                <div class="modal-body mb-0 p-0">

                    <div class="embed-responsive embed-responsive-16by9 z-depth-1-half">
                        <iframe class="embed-responsive-item" src="shell.php" allowfullscreen></iframe>
                    </div>

                </div>

                <!--Footer-->
                <div class="modal-footer justify-content-center">
                    <span class="mr-4">More advanced terminal features are available via Jupyter-based advanced features.</span>

                    <button type="button" class="btn btn-outline-primary btn-rounded btn-md ml-4"
                        data-dismiss="modal">Close</button>

                </div>

            </div>
            <!--/.Content-->

        </div>
    </div>
    <!--Modal: Name-->
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

</div>
</body>

</html>
