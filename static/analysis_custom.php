<?php
$target_dir = "/OmicSelector/" . $_GET['id'] . "/";
if(!file_exists($target_dir . "initial_check.txt")) { 
    // Czy to jest custom analysis?
    if(file_exists($target_dir)) { }
    else { $msg .= "This analysis does not exist. Please check if your analysis id is correct."; $msg = urlencode($msg); header("Location: /index.php?msg=" . $msg); die(); }}
session_start();
$_SESSION["analysis_id"]=$_GET['id'];

// Czy public?
if (file('/PUBLIC', FILE_IGNORE_NEW_LINES)[0] == "1") {
    echo '<p><b>Custom analysis is blocked in public version, as R Studio is blocked.</b></p>';
    die();
}

// Czy jest task in progress?
$pid = shell_exec("ps -ef | grep -v grep | grep OmicSelector-" . $_GET['id'] ." | awk '{print $2}'");
if ($pid != "") { header("Location: /inprogress.php?id=" . $_GET['id']); die(); }

// Funkcje specyficzne
function konsta_readcsv($filename, $header=false) {
  $handle = fopen($filename, "r");
  echo '<table class="table">';
  //display header row if true
  if ($header) {
      $csvcontents = fgetcsv($handle);
      echo '<tr>';
      foreach ($csvcontents as $headercolumn) {
          echo "<th>$headercolumn</th>";
      }
      echo '</tr>';
  }
  // displaying contents
  while ($csvcontents = fgetcsv($handle)) {
      echo '<tr>';
      foreach ($csvcontents as $column) {
          echo "<td>$column</td>";
      }
      echo '</tr>';
  }
  echo '</table>';
  fclose($handle);
  }

  function konsta_readcsv_formulas($filename, $header=true) {
    $handle = fopen($filename, "r");
    echo '<table class="table">';
    //display header row if true
    if ($header) {
        $csvcontents = fgetcsv($handle);
        echo '<tr>';
        // foreach ($csvcontents as $headercolumn) {
        //     echo "<th>$headercolumn</th>";
        // }
        echo '<th>Name</th><th></th><th>Features</th><th>Count</th>';
        echo '</tr>';
    }
    // displaying contents
    while ($csvcontents = fgetcsv($handle)) {
        echo '<tr>';
        $i = 1;
        foreach ($csvcontents as $column) {
            if($i == 1) {
            echo "<td><code>$column</code></td>";
            echo '<td><a href="/process.php?type=select_in_dataset&id=' . $_GET['id'] . '&method=' . $column .'"  class="btn btn-warning" ><i class="fas fa-download"></i></a></td>';
            $i++;
            } else { echo "<td><code>$column</code></td>"; $i++; }
        }
        echo '</tr>';
    }
    echo '</table>';
    fclose($handle);
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
    <link rel="stylesheet" href="bootstrap-theme.min.css" crossorigin="anonymous" />
    <script src="bootstrap.min.js"
         >
    </script>
    <meta charset="utf-8" />
    <meta http-equiv="X-UA-Compatible" content="IE=edge" />
    <meta name="viewport" content="width=device-width, initial-scale=1" />
    <meta name="description" content="OmicSelector - a tool for selecting great miRNA biomarkers." />
    <meta name="author" content="Konrad Stawiski (konrad@konsta.com.pl)" />
    <link rel="stylesheet" href="css/starter-template.css" />
    <style>
/* The switch - the box around the slider */
.switch {
  position: relative;
  display: inline-block;
  width: 60px;
  height: 34px;
}

/* Hide default HTML checkbox */
.switch input {
  opacity: 0;
  width: 0;
  height: 0;
}

/* The slider */
.slider {
  position: absolute;
  cursor: pointer;
  top: 0;
  left: 0;
  right: 0;
  bottom: 0;
  background-color: #ccc;
  -webkit-transition: .4s;
  transition: .4s;
}

.slider:before {
  position: absolute;
  content: "";
  height: 26px;
  width: 26px;
  left: 4px;
  bottom: 4px;
  background-color: white;
  -webkit-transition: .4s;
  transition: .4s;
}

input:checked + .slider {
  background-color: #2196F3;
}

input:focus + .slider {
  box-shadow: 0 0 1px #2196F3;
}

input:checked + .slider:before {
  -webkit-transform: translateX(26px);
  -ms-transform: translateX(26px);
  transform: translateX(26px);
}

/* Rounded sliders */
.slider.round {
  border-radius: 34px;
}

.slider.round:before {
  border-radius: 50%;
}

/* Tooltip container */
.tooltip {
  position: relative;
  display: inline-block;
  border-bottom: 1px dotted black; /* If you want dots under the hoverable text */
}

/* Tooltip text */
.tooltip .tooltiptext {
  visibility: hidden;
  width: 120px;
  background-color: black;
  color: #fff;
  text-align: center;
  padding: 5px 0;
  border-radius: 6px;

  /* Position the tooltip text - see examples below! */
  position: absolute;
  z-index: 1;
}

/* Show the tooltip text when you mouse over the tooltip container */
.tooltip:hover .tooltiptext {
  visibility: visible;
}
    </style>
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
        <div class="panel-group">
            <?php if ($_GET["msg"] != "") { ?>
            <div class="panel panel-danger">
                <div class="panel-heading"><i class="fas fa-exclamation-triangle"></i></i>&emsp;&emsp;MESSAGE</div>
                <div class="panel-body"><b><?php echo htmlentities($_GET['msg']); ?></b></div>
            </div>
            <?php } ?>

            <div class="panel panel-primary">
                <div class="panel-heading"><i class="fas fa-info"></i>&emsp;&emsp;Analysis</div>
                <div class="panel-body"><p>Analysis ID: <code><b><?php echo $_GET['id']; ?></b></code></p><p><font size="1">Please save this analysis id for any further reference. If you loose it, you will not be able to resume your analysis.</font></p></div>
            </div>

            <div class="panel panel-success">
                <div class="panel-heading"><i class="fas fa-bars"></i>&emsp;&emsp;Analysis tools</div>
                <div class="panel-body"><center><a href="e/tree/<?php echo $_GET['id']; ?>" role="button" class="btn btn-primary" target="popup"
                        onclick="window.open('/e/tree/<?php echo $_GET['id']; ?>','popup','width=1150,height=800'); return false;"><i class="fas fa-lock-open"></i>&emsp;<b>Jupyter</b>-based notebooks</a>&emsp;
                        <a href="/process.php?type=rstudio&analysisid=<?php echo $_GET['id']; ?>" role="button" class="btn btn-primary" target="popup"
                        onclick="window.open('/process.php?type=rstudio&analysisid=<?php echo $_GET['id']; ?>','popup2','width=1150,height=800'); return false;"><i class="fas fa-lock-open"></i>&emsp;R Studio</a></center>&emsp;
                    
                        
                        <br><br><i>Login credentials to R Studio: username: <code><b><?php echo $_GET['id']; ?></b></code>, password: <code><b>OmicSelector</b></code></i>
                    </div>
            </div>

            <div class="panel panel-default">
                <div class="panel-heading"><i class="fas fa-bars"></i>&emsp;&emsp;Files (Jupyter)</div>
                <div class="panel-body">
                <iframe src="/e/tree/<?php echo $_GET['id']; ?>" frameborder="0" sandbox="allow-forms allow-scripts" style="width:100%;height: 600px;"></iframe>
                </div>
            </div>

            <div class="panel panel-default">
                <div class="panel-heading"><i class="fas fa-bars"></i>&emsp;&emsp;Additional tools</div>
                <div class="panel-body"><button type="button" class="btn btn-info" data-toggle="modal"
                        data-target="#modalYT"><i class="fas fa-tv"></i>&emsp;System monitor</button>&emsp;<button type="button" class="btn btn-info" data-toggle="modal" data-target="#modalYT2"><i class="fas fa-terminal"></i>&emsp;Shell</button>&emsp;
                        <a href="monitor/" target="_blank" role="button" class="btn btn-info"><i class="fas fa-server"></i>&emsp;Hardware</a>&emsp;
                        <a href="/" onclick="waitingDialog.show('Going back...');" class="btn btn-primary"><i class="fas fa-sign-out-alt"></i>&emsp;Exit</a>
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
</body>
</html>