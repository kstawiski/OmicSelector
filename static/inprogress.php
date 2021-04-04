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
        $('#log').animate({scrollTop:document.getElementById("log").scrollHeight}, 'slow');
    }
});</script>
	<?php } ?>
    </div></div>
    <div class="panel panel-default">
                <div class="panel-heading"><i class="fas fa-file-medical-alt"></i>&emsp;&emsp;Task details</div>
                <div class="panel-body">
    <p><b>Task progress (log file, last 1000 lines):</b> (this is updated in the real-time)</p>
    <p></p><pre id="log"></pre></p>
    <p><a href="/e/files/<?php echo $_GET['id']; ?>/task.log" class="btn btn-primary"><i class="fas fa-file-download"></i>&emsp;Download whole log file (text file)</a>&emsp;
<?php if($skonczone == 1) { ?><a href="/analysis.php?id=<?php echo $analysis_id; ?>" onclick="waitingDialog.show('Loading...');" class="btn btn-success"><i class="fas fa-chart-line"></i>&emsp;Go forward</a><?php } ?>
    </p>
    <hr>
    <p>Process details:</p>
    <p><pre><?php echo $task_process; ?></pre></p>
    </div></div>

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


    <div class="panel panel-default">
                <div class="panel-heading"><i class="fas fa-bars"></i>&emsp;&emsp;Additional tools</div>
                <div class="panel-body"><button type="button" class="btn btn-info" data-toggle="modal"
                        data-target="#modalYT"><i class="fas fa-tv"></i>&emsp;CPU Monitor</button>&emsp;
                        <button type="button" class="btn btn-info" data-toggle="modal"
                        data-target="#modalYT2"><i class="fas fa-microchip"></i>&emsp;GPU Monitor</button>&emsp;
                        <a href="e/tree/<?php echo $_GET['id']; ?>" role="button" class="btn btn-primary" target="popup"
                        onclick="window.open('/e/tree/<?php echo $_GET['id']; ?>','popup','width=1150,height=800'); return false;"><i class="fas fa-lock-open"></i>&emsp;Advanced features (Jupyter)</a>&emsp;
                        <a href="/process.php?type=rstudio&analysisid=<?php echo $_GET['id']; ?>" role="button" class="btn btn-primary" target="popup"
                        onclick="window.open('/process.php?type=rstudio&analysisid=<?php echo $_GET['id']; ?>','popup','width=1150,height=800'); return false;"><i class="fas fa-lock-open"></i>&emsp;Advanced features (R Studio)</a>&emsp;
                        <a href="/process.php?type=radiant&analysisid=<?php echo $_GET['id']; ?>" role="button" class="btn btn-primary" target="popup"
                        onclick="window.open('/process.php?type=radiant&analysisid=<?php echo $_GET['id']; ?>','popup3','width=1150,height=800'); return false;"><i class="fas fa-lock-open"></i>&emsp;Radiant</a>&emsp;
                        <a href="/process.php?type=vscode&analysisid=<?php echo $_GET['id']; ?>" role="button" class="btn btn-primary" target="popup"
                        onclick="window.open('/process.php?type=vscode&analysisid=<?php echo $_GET['id']; ?>','popup3','width=1150,height=800'); return false;"><i class="fas fa-lock-open"></i>&emsp;VS Code</a>
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
                        <iframe class="embed-responsive-item" src="gpu.php" allowfullscreen></iframe>
                    </div>

                </div>

                <!--Footer-->
                <div class="modal-footer justify-content-center">
                    <span class="mr-4">Monitoring GPU via <code>nvidia-smi</code> every 5 seconds.</span>

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
  </body>
</html>