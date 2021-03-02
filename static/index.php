<?php
exec("chown -R OmicSelector /OmicSelector");
require_once 'class.formr.php';
if (!file_exists('/OmicSelector/var_status.txt')) { file_put_contents('/OmicSelector/var_status.txt', "[0] INITIAL (UNCONFIGURED)"); } // Wyjściowy status.
$status = file_get_contents('/OmicSelector/var_status.txt');
if (file_exists('/version.txt')) { $version = file_get_contents('/version.txt'); } else { $version = "init"; }
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
    <meta name="author" content="Konrad Stawiski (konrad.stawiski@umed.lodz.pl)" />
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

$(document).on('click', '.panel div.clickable', function (e) {
    var $this = $(this); //Heading
    var $panel = $this.parent('.panel');
    var $panel_body = $panel.children('.panel-body');
    var $display = $panel_body.css('display');

    if ($display == 'block') {
        $panel_body.slideUp();
    } else if($display == 'none') {
        $panel_body.slideDown();
    }
});

$(document).ready(function(e){
    var $classy = '.panel.autocollapse';

    var $found = $($classy);
    $found.find('.panel-body').hide();
    $found.removeClass($classy);
});

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
    <style>.clickable { cursor: pointer; }</style>
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
                <div class="panel-body"><p>Welcome to <b>OmicSelector</b> - the software intended to find the best biomarker signature based on NGS or other omic modalities (miRNA-seq, RNA-seq).</p>

        <p>Your current version of software: <code>OmicSelector v1.0.<?php echo $version; ?></code></p>
        <?php
            if (file('/PUBLIC', FILE_IGNORE_NEW_LINES)[0] == "1") {
                echo '<p><b>This is a public (demo) version of software. Some hacking options are restricted. The public version cleans data and restarts every Monday (at night).</b></p>';
            }
        ?>
        </div>
            </div>
            <div class="panel panel-primary">
                <div class="panel-heading"><i class="fas fa-chart-pie"></i>&emsp;&emsp;Analysis</div>
                <div class="panel-body">
                <p>You can start new analysis or resume your previous one.</p>
                <table class="table">
                    <tr><td width="50%">
                    <h4>Start new analysis:</h4>
                    <p><a href="/start.php" role="button" class="btn btn-primary"><i class="fas fa-plus"></i>&emsp;New OmicSelector analysis</a></p>
                    <p><a href="/start_custom.php" role="button" class="btn btn-warning"><i class="fas fa-plus"></i>&emsp;New custom analysis</a></p>
                    </td>
                    <td width="50%">
                        <h4>Resume analysis:</h4>
                <p><form action="/analysis.php" method="get">
                    <p>
                    <?php
                    if (file('/PUBLIC', FILE_IGNORE_NEW_LINES)[0] != "1") {
                    ?>
                    <link href="https://cdn.jsdelivr.net/npm/select2@4.1.0-beta.1/dist/css/select2.min.css" rel="stylesheet" />
                    <script src="https://cdn.jsdelivr.net/npm/select2@4.1.0-beta.1/dist/js/select2.min.js"></script>
                    <script>$(document).ready(function () { $('#id').select2(); }); </script>
                    <select id="id" name="id" placeholder="Provide analysis ID" class="form-control">
                        <?php
                        $directories = glob('/OmicSelector' . '/*' , GLOB_ONLYDIR);

                        foreach($directories as $item){
                            $item2 = basename($item);
                            if($item2 != "OmicSelector") {
                            echo "<option value='$item2'>$item2</option>"; }
                        }
                        ?>
                    </select>
                    <?php } else { ?>
                    <input type="text" id="id" name="id" placeholder="Provide analysis ID" class="form-control" value="<?php echo $prev_analysis; ?>" autocomplete="analysisid name">
                    <?php } ?>
                    
                    
                    </p>
                <p><button type="submit" class="btn btn-success" value="Upload" name="submit" onclick="waitingDialog.show('Loading...');"><i class="fas fa-folder-open"></i>&emsp;Resume analysis</button></p>
                </form></p>
            </td></tr>
            </table>
                </div>
            </div>


            <div class="panel panel-warning">
                <div class="panel-heading"><i class="fas fa-puzzle-piece"></i>&emsp;&emsp;Preprocessing and postprocessing extensions</div>
                <div class="panel-body">

                <div class="panel panel-default autocollapse">
                        <div class="panel-heading clickable">
                            <h3 class="panel-title" id="deep_learning">
                            <i class="fas fa-code-branch"></i>&emsp;Differential expression analysis using corrected t-test.
                            </h3>
                        </div>
                        <div class="panel-body">
                            <p>This app allows you to conduct DE analysis. The file should be prepared in the same as for the OmicSelector analysis.</p>
                            <a href="/tool_de/" role="button" class="btn btn-primary" target="popup" onclick="window.open('/tool_de/','popup3','width=1150,height=800'); return false;"><i class="fas fa-external-link-alt"></i>&emsp;Open tool</a>&emsp;
                        </div>
                    </div>
                
                    <div class="panel panel-default autocollapse">
                        <div class="panel-heading clickable">
                            <h3 class="panel-title" id="deep_learning">
                            <i class="fas fa-code-branch"></i>&emsp;Imputation of missing data.
                            </h3>
                        </div>
                        <div class="panel-body">
                            <p>This app allows you to impute missing data using predictive mean matching or mean. The file should be prepared in the same as for the OmicSelector analysis.</p>
                            <a href="/tool_impute/" role="button" class="btn btn-primary" target="popup" onclick="window.open('/tool_impute/','popup3','width=1150,height=800'); return false;"><i class="fas fa-external-link-alt"></i>&emsp;Open tool</a>&emsp;
                        </div>
                    </div>

                
                    <div class="panel panel-default autocollapse">
                        <div class="panel-heading clickable">
                            <h3 class="panel-title" id="deep_learning">
                            <i class="fas fa-code-branch"></i>&emsp;Correct batch effect using ComBat.
                            </h3>
                        </div>
                        <div class="panel-body">
                            <p>This app allows you to correct batch effect using ComBat. The file should be prepared in the same as for the OmicSelector analysis.</p>
                            <a href="/tool_batch/" role="button" class="btn btn-primary" target="popup" onclick="window.open('/tool_batch/','popup3','width=1150,height=800'); return false;"><i class="fas fa-external-link-alt"></i>&emsp;Open tool</a>&emsp;
                        </div>
                    </div>

                    <div class="panel panel-default autocollapse">
                        <div class="panel-heading clickable">
                            <h3 class="panel-title" id="deep_learning">
                            <i class="fas fa-code-branch"></i>&emsp;Generate heatmap for exploratory analysis.
                            </h3>
                        </div>
                        <div class="panel-body">
                            <p>The file should be prepared in the same as for the OmicSelector analysis.</p>
                            <a href="/tool_heatmap/" role="button" class="btn btn-primary" target="popup" onclick="window.open('/tool_heatmap/','popup3','width=1150,height=800'); return false;"><i class="fas fa-external-link-alt"></i>&emsp;Open tool</a>&emsp;
                            
                        </div>
                    </div>


                    <div class="panel panel-default autocollapse">
                        <div class="panel-heading clickable">
                            <h3 class="panel-title" id="deep_learning">
                            <i class="fas fa-code-branch"></i>&emsp;<b>[DEEP LEARNING]</b> Predict with developed deep learning model.
                            </h3>
                        </div>
                        <div class="panel-body">
                            <p>This step requires the model zip file.</p>
                            <a href="/deeplearning_model/" role="button" class="btn btn-primary" target="popup" onclick="window.open('/deeplearning_model/','popup3','width=1150,height=800'); return false;"><i class="fas fa-external-link-alt"></i>&emsp;Open tool</a>&emsp;
                            
                        </div>
                    </div>
                </div>
            </div>



        <div class="panel panel-default">
                <div class="panel-heading"><i class="fas fa-bars"></i>&emsp;&emsp;Additional tools</div>
                <div class="panel-body"><button type="button" class="btn btn-info" data-toggle="modal"
                        data-target="#modalYT"><i class="fas fa-tv"></i>&emsp;System monitor</button>&emsp;<button type="button" class="btn btn-info" data-toggle="modal" data-target="#modalYT2"><i class="fas fa-terminal"></i>&emsp;Shell</button>&emsp;
                        <a href="process.php?type=init_update" role="button" onclick="waitingDialog.show('Starting update...');" class="btn btn-primary"><i class="fas fa-arrow-up"></i></i>&emsp;Update</a>
                    &emsp;<a href="e/notebooks/OmicSelector/vignettes/Tutorial.Rmd" role="button" class="btn btn-primary" target="_blank"><i class="fas fa-graduation-cap"></i>&emsp;Learn R package</a>
                    <?php
                    // Advanced features in private version:
                    if (file('/PUBLIC', FILE_IGNORE_NEW_LINES)[0] != "1") {
                        echo '&emsp;<a href="e/tree" role="button" class="btn btn-primary" target="_blank"><i class="fas fa-puzzle-piece"></i>&emsp;Advanced features</a>';
                    }
                    ?>
                    </div>
        </div>


        <script type="text/javascript" src="monitor/gauge/jquery-asPieProgress.js"></script>
        <script type="text/javascript">
            $(document).ready(function () {
                // Example with grater loading time - loads longer
                $('.pie_progress_temperature,.pie_progress_cpu, .pie_progress_mem, .pie_progress_disk').asPieProgress({});
		getTemp();
                getCpu();
                getMem();
                getDisk();
            });

            function getTemp() {
                $.ajax({
                    url: 'monitor/temperature.json.php',
                    success: function (response) {
                        update('temperature', response);
                        setTimeout(function () {
                            getTemp();
                        }, 1000);
                    }
                });
            }


            function getCpu() {
                $.ajax({
                    url: 'monitor/cpu.json.php',
                    success: function (response) {
                        update('cpu', response);
                        setTimeout(function () {
                            getCpu();
                        }, 1000);
                    }
                });
            }

            function getMem() {
                $.ajax({
                    url: 'monitor/memory.json.php',
                    success: function (response) {
                        update('mem', response);

                        setTimeout(function () {
                            getMem();
                        }, 1000);
                    }
                });
            }

            function getDisk() {
                $.ajax({
                    url: 'monitor/disk.json.php',
                    success: function (response) {
                        update('disk', response);
                        setTimeout(function () {
                            getDisk();
                        }, 1000);
                    }
                });
            }

            function update(name, response) {
                $('.pie_progress_' + name).asPieProgress('go', response.percent);
                $("#" + name + "Div div.title").text(response.title);
                //$("#" + name + "Div pre").text(response.output.join('\n'));
            }
        </script>
        <link rel="stylesheet" href="monitor/gauge/css/asPieProgress.css">
        <p></p>
        <div class="panel panel-default">
                <div class="panel-heading"><i class="fas fa-heartbeat"></i>&emsp;&emsp;Resources (monitor)</div>
                <div class="panel-body">
                <div class="col-xs-3 col-sm-3 col-lg-3" id="cpuDiv">                        
                <div class="pie_progress_cpu" role="progressbar" data-goal="33">
                    <div class="pie_progress__number">0%</div>
                    <div class="pie_progress__label">CPU</div>
                </div>
  
                <div class='title'></div>
            </div>
            <div class="col-xs-3 col-sm-3 col-lg-3" id="memDiv">
                <div class="pie_progress_mem" role="progressbar" data-goal="33">
                    <div class="pie_progress__number">0%</div>
                    <div class="pie_progress__label">Memory</div>
                </div>
  
                <div class='title'></div>
            </div>
            <div class="col-xs-3 col-sm-3 col-lg-3" id="diskDiv">
                <div class="pie_progress_disk" role="progressbar" data-goal="33">
                    <div class="pie_progress__number">0%</div>
                    <div class="pie_progress__label">Disk</div>
                </div>

                <div class='title'></div>
            </div>
            <div class="col-xs-3 col-sm-3 col-lg-3" id="temperatureDiv">
                <div class="pie_progress_temperature" role="progressbar" data-goal="33">
                    <div class="pie_progress__number">0°</div>
                    <div class="pie_progress__label">Temperature</div>
                </div>
                <div class='title' style="display:none;"></div>
            </div>



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
                    class="fas fa-envelope"></i> konrad.stawiski@umed.lodz.pl&emsp;&emsp;&emsp;<i
                    class="fas fa-globe-europe"></i>
                 <a href="https://biostat.umed.pl/OmicSelector/"
                    target="_blank">https://biostat.umed.pl/OmicSelector/</a></span>
                    <p>&emsp;</p>
        </div>
    </footer>
    <!-- /.container -->

</div>
</body>

</html>
