<html>
<?php 
if(!isset($_GET['id'])) {
$czy_dziala = "Not running";
$target_log = "/task.log";
if(file_exists($target_log)) { 
    $pid = shell_exec("ps -ef | grep -v grep | grep OmicSelector-task | awk '{print $2}'");
    $zawartosc_logu = file_get_contents($target_log);
    $task_process = shell_exec('ps -ef | grep -v grep | grep OmicSelector-task');
    if($task_process == "") { $task_process = "NOT RUNNING"; }
    // $czy_dziala = shell_exec('ps -ef | grep -v grep | grep OmicSelector-task | wc -l');
    if ($pid != "") { $czy_dziala = "Running"; }
    $skonczone = 0; if (strpos($zawartosc_logu, '[OmicSelector: TASK COMPLETED]') !== false) { $skonczone = 1; } 
    // if (strpos($zawartosc_logu, 'Error') !== false) { $skonczone = 1; } 
} else { $msg = urlencode("The task was not initialized. Please run it again."); header("Location: /?msg=" . $msg); die(); }
}
else {
    $czy_dziala = "Not running";
    $analysis_id = $_GET['id'];
    $target_dir = "/OmicSelector/" . $analysis_id . "/";
    $target_log = $target_dir . "task.log";
    if(file_exists($target_log)) { 
        $pid = shell_exec("ps -ef | grep -v grep | grep OmicSelector-" . $analysis_id ." | awk '{print $2}'");
        $zawartosc_logu = file_get_contents($target_log);
        $task_process = shell_exec('ps -ef | grep -v grep | grep OmicSelector-' . $analysis_id);
        if($task_process == "") { $task_process = "NOT RUNNING"; }
        // $czy_dziala = shell_exec('ps -ef | grep -v grep | grep OmicSelector-task | wc -l');
        if ($pid != "") { $czy_dziala = "Running"; }
        $skonczone = 0; if (strpos($zawartosc_logu, '[OmicSelector: TASK COMPLETED]') !== false) { $skonczone = 1; } 
        // if (strpos($zawartosc_logu, 'Error') !== false) { $skonczone = 1; } 
} else { $msg = urlencode("The task was not initialized. Please run it again. Go again to your analysis and init it."); header("Location: /?msg=" . $msg); die(); }
}
// if()
?>

<head>
    <title>OmicSelector (task in progress)</title>
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
	<script src="konsta.js"></script>
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
            <div class="panel panel-primary">
                <div class="panel-heading"><i class="fas fa-info"></i>&emsp;&emsp;Task status</div>
                <div class="panel-body">
    <?php if($skonczone == 1) { ?>
	<p id="msg"><b>The task is finished.</b> Please go to the next step to analyze the results and move to next steps.</p>
    <p><a href="/analysis.php?id=<?php echo $analysis_id; ?>" onclick="waitingDialog.show('Loading...');" class="btn btn-success"><i class="fas fa-chart-line"></i>&emsp;Go forward</a></p>
    <script>
        $( document ).ready(function() {
        $(document).scrollTop($(document).height()); 
        get_log(); 
    
    function get_log(){
        var feedback = $.ajax({
            type: "GET",
            url: "api.php",
            async: false,
            data: {
                "typ" : "odczyt",
                "plik" : "<?php echo $target_log; ?>"
            },
            success: function() {
                // setTimeout(function(){ get_log();}, 1000);
            }
        }).responseText;
    $('#log').html(feedback);
    }
});
    </script>
	<?php } else { ?>
		<p id="msg"><b>Status:</b><pre><?php echo $czy_dziala . " PID=" . $pid; ?></pre>Please wait and do not use the app until this is finished. </p>
        <?php if($skonczone == 0 && $pid == "") {
            echo "<p style=\"color:red;\"><b><i class=\"fas fa-exclamation-triangle\"></i> The task stopped. Probably due to error. Please check out the log below, cancel the task using the button below and try to fix the issue. If you think this is a bug please report it on GitHub.</b></p>";
        } ?>
        <p><a href="/process.php?type=cancel&pid=<?php echo $pid; ?>" onclick="return confirm('Are you sure? This will kill the process of task processing and you will need to reconfigure this task!')" class="btn btn-danger"><i class="fas fa-skull-crossbones"></i>&emsp;Cancel task</a></p>
        <script>
        $( document ).ready(function() {
    $(document).scrollTop($(document).height()); 
    setInterval(get_log, 2000);
    
    
    function get_log(){
        var feedback = $.ajax({
            type: "GET",
            url: "api.php",
            async: false,
            data: {
                "typ" : "odczyt",
                "plik" : "<?php echo $target_log; ?>"
            },
            success: function() {
                //setTimeout(function(){ get_log();}, 1000);
            }
        }).responseText;
        
        var substring = "[OmicSelector: TASK COMPLETED]";
        console.log("Refreshing log..");
        if(feedback.includes(substring)) { location.reload(); }

        $('#log').html(feedback);
    }
});</script>
	<?php } ?>
    </div></div>
    <div class="panel panel-default">
                <div class="panel-heading"><i class="fas fa-file-medical-alt"></i>&emsp;&emsp;Task details</div>
                <div class="panel-body">
    <p><b>Task progress (log file):</b> (this is updated in the real-time)</p>
    <p></p><pre id="log"></pre></p>
    <hr>
    <p>Process details:</p>
    <p><pre><?php echo $task_process; ?></pre></p>
    </div></div>
    <div class="panel panel-default">
                <div class="panel-heading"><i class="fas fa-bars"></i>&emsp;&emsp;Additional tools</div>
                <div class="panel-body"><button type="button" class="btn btn-info" data-toggle="modal"
                        data-target="#modalYT"><i class="fas fa-tv"></i>&emsp;System monitor</button>&emsp;<button type="button" class="btn btn-info" data-toggle="modal" data-target="#modalYT2"><i class="fas fa-terminal"></i>&emsp;Shell</button>&emsp;
                        <a href="monitor/" target="_blank" role="button" class="btn btn-info"><i class="fas fa-server"></i>&emsp;Hardware</a>&emsp;<a href="e/tree/<?php echo $_GET['id']; ?>" role="button" class="btn btn-primary" target="popup"
                        onclick="window.open('/e/tree/<?php echo $_GET['id']; ?>','popup','width=1150,height=800'); return false;"><i class="fas fa-lock-open"></i>&emsp;Advanced features (Jupyter)</a>&emsp;
                        <a href="/process.php?type=rstudio&analysisid=<?php echo $_GET['id']; ?>" role="button" class="btn btn-primary" target="popup"
                        onclick="window.open('/process.php?type=rstudio&analysisid=<?php echo $_GET['id']; ?>','popup','width=1150,height=800'); return false;"><i class="fas fa-lock-open"></i>&emsp;Advanced features (R Studio)</a>&emsp;
                        <a href="/process.php?type=radiant&analysisid=<?php echo $_GET['id']; ?>" role="button" class="btn btn-primary" target="popup"
                        onclick="window.open('/process.php?type=radiant&analysisid=<?php echo $_GET['id']; ?>','popup3','width=1150,height=800'); return false;">Radiant</a>&emsp;
                        <a href="/process.php?type=vscode&analysisid=<?php echo $_GET['id']; ?>" role="button" class="btn btn-primary" target="popup"
                        onclick="window.open('/process.php?type=vscode&analysisid=<?php echo $_GET['id']; ?>','popup3','width=1150,height=800'); return false;">Radiant</a>
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
                    class="fas fa-envelope"></i> konrad@konsta.com.pl&emsp;&emsp;&emsp;<i
                    class="fas fa-globe-europe"></i>
                <a href="https://biostat.umed.pl" taret="_blank">https://biostat.umed.pl</a>&emsp;&emsp;&emsp;<i
                    class="fab fa-github"></i> <a href="https://kstawiski.github.io/OmicSelector/"
                    target="_blank">https://kstawiski.github.io/OmicSelector/</a></span>
                    <p>&emsp;</p>
        </div>
    </footer>
  </body>
</html>