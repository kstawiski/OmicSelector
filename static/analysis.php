<?php
$target_dir = "/OmicSelector/" . $_GET['id'] . "/";
if (!file_exists($target_dir . "initial_check.txt")) {
    // Czy to jest custom analysis?
    if (file_exists($target_dir)) {
        header("Location: /analysis_custom.php?id=" . $_GET['id']);
        die();
    } else {
        $msg .= "This analysis does not exist. Please check if your analysis id is correct.";
        $msg = urlencode($msg);
        header("Location: /index.php?msg=" . $msg);
        die();
    }
}
session_start();
$_SESSION["analysis_id"] = $_GET['id'];


// Czy jest task in progress?
if($_GET['override'] != 1) {
$pid = shell_exec("ps -ef | grep -v grep | grep OmicSelector-" . $_GET['id'] . " | awk '{print $2}'");
if ($pid != "") {
    header("Location: /inprogress.php?id=" . $_GET['id']);
    die();
}
}

// Funkcje specyficzne
function konsta_readcsv($filename, $header = false)
{
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

function konsta_readcsv_formulas($filename, $header = true)
{
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
            if ($i == 1) {
                echo "<td><code>$column</code></td>";
                echo '<td style="white-space: nowrap !important;"><a href="/process.php?type=select_in_dataset&id=' . $_GET['id'] . '&method=' . $column . '"  class="btn btn-warning" ><i class="fas fa-download"></i></a>&nbsp;<a href="/process.php?type=delete_set_from_final&id=' . $_GET['id'] . '&name=' . $column . '"  class="btn btn-danger" onclick="return confirm(\'Are you sure that you want to delete this feature set from benchmarking?\')"><i class="fas fa-trash"></i></a></td>';
                $i++;
            } else {
                echo "<td><code>$column</code></td>";
                $i++;
            }
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
    <script src="https://code.jquery.com/jquery-3.5.1.js" crossorigin="anonymous"></script>
    <script src="jquery-ui.js" type="text/javascript"></script>
    <link rel="stylesheet" href="bootstrap.min.css" />
    <link rel="stylesheet" href="bootstrap-theme.min.css" crossorigin="anonymous" />
    <script src="bootstrap.min.js">
    </script>
    <meta charset="utf-8" />
    <meta http-equiv="X-UA-Compatible" content="IE=edge" />
    <meta name="viewport" content="width=device-width, initial-scale=1" />
    <meta name="description" content="OmicSelector - a tool for selecting great miRNA biomarkers." />
    <meta name="author" content="Konrad Stawiski (konrad.stawiski@umed.lodz.pl)" />
    <link rel="stylesheet" href="css/starter-template.css" />
    <style>.clickable { cursor: pointer; }</style>
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

        input:checked+.slider {
            background-color: #2196F3;
        }

        input:focus+.slider {
            box-shadow: 0 0 1px #2196F3;
        }

        input:checked+.slider:before {
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
            border-bottom: 1px dotted black;
            /* If you want dots under the hoverable text */
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
    <script src="https://cdnjs.cloudflare.com/ajax/libs/font-awesome/5.14.0/js/all.min.js"></script>
    <!-- Global site tag (gtag.js) - Google Analytics -->
    <script async src="https://www.googletagmanager.com/gtag/js?id=UA-53584749-8"></script>
    <script>
        window.dataLayer = window.dataLayer || [];

        function gtag() {
            dataLayer.push(arguments);
        }
        gtag('js', new Date());

        gtag('config', 'UA-53584749-8');
    </script>
    <script type="text/javascript">
        $(".btn-success").click(function(event) {
            waitingDialog.show('Processing.. Please wait...');
        });

        $(document).on('click', '.panel div.clickable', function(e) {
            var $this = $(this); //Heading
            var $panel = $this.parent('.panel');
            var $panel_body = $panel.children('.panel-body');
            var $display = $panel_body.css('display');

            if ($display == 'block') {
                $panel_body.slideUp();
            } else if ($display == 'none') {
                $panel_body.slideDown();
            }
        });

        $(document).ready(function(e) {
            var $classy = '.panel.autocollapse';

            var $found = $($classy);
            $found.find('.panel-body').hide();
            $found.removeClass($classy);
        });

        var waitingDialog = waitingDialog || (function($) {
            'use strict';

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
                show: function(message, options) {
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
                        $dialog.off('hidden.bs.modal').on('hidden.bs.modal', function(e) {
                            settings.onHide.call($dialog);
                        });
                    }
                    // Opening dialog
                    $dialog.modal();
                },
                /**
                 * Closes dialog
                 */
                hide: function() {
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
                <div class="panel-body">
                    <p>Analysis ID: <code><b><?php echo $_GET['id']; ?></b></code></p>
                    <p>
                        <font size="1">Please save this analysis id for any further reference. If you loose it, you will not be able to resume your analysis.</font>
                    </p>
                </div>
            </div>

            <div class="panel panel-info">
                <div class="panel-heading"><i class="fas fa-table"></i>&emsp;&emsp;Data</div>
                <div class="panel-body">
                    <p>Initial check status:
                        <pre style="white-space: pre-wrap;"><?php echo file_get_contents($target_dir . "initial_check.txt"); ?></pre>
                    </p>
                    <p>Initial check status:
                        <code><b><?php $var_initcheck = file_get_contents($target_dir . 'var_initcheck.txt');
                                    echo $var_initcheck; ?></b></code>
                    </p>
                    <p><a href="viewer.php?f=<?php echo $_GET['id']; ?>/data_start.csv" class="btn btn-info" role="button" target="popup" onclick="window.open('viewer.php?f=<?php echo $_GET['id']; ?>/data_start.csv','popup','width=1150,height=800'); return false;">View
                            data</a></p>
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

            <div class="panel panel-primary">
                <div class="panel-heading"><i class="fas fa-chart-bar"></i>&emsp;&emsp;Exploratory analysis</div>
                <div class="panel-body">
                    <?php
                    $images = glob($target_dir . "exploratory_*.png");
                    foreach ($images as $image) {
                        $image3 = str_replace("/OmicSelector", "/e/files", $image);
                        $image2 = str_replace("/OmicSelector", "/e/view", $image);
                        echo '<a href="' . $image2 . '" target="_blank" onclick="window.open(\'' . $image2 . '\',\'popup\',\'width=600,height=600\'); return false;"><img src="' . $image3 . '" width="49%" /></a>';
                    }
                    ?>
                    <p>
                        <font size="1">Notes: <i>If mix is not labeled the heatmap was constructed based on training set. Some of the heatmaps use raw expression, some using z-scoring. Features marked on vulcano plot are significant in DE. You can re-create and customize those plots below.</i></font>
                    </p>
                </div>
                <table class="table">

                    <tbody>

                        <tr>
                            <td>Training set:</td>
                            <td><a href="viewer.php?f=<?php echo $_GET['id']; ?>/mixed_train.csv" class="btn btn-info" role="button" target="popup" onclick="window.open('viewer.php?f=<?php echo $_GET['id']; ?>/mixed_train.csv','popup','width=1150,height=800'); return false;"><i class="fas fa-search-plus"></i> View</a>&emsp;<a href="/e/files/<?php echo $_GET['id']; ?>/mixed_train.csv" class="btn btn-warning"><i class="fas fa-download"></i> Download</a></td>
                        </tr>

                        <tr>
                            <td>Training set (balanced):</td>
                            <td><a href="viewer.php?f=<?php echo $_GET['id']; ?>/mixed_train.csv" class="btn btn-info" role="button" target="popup" onclick="window.open('viewer.php?f=<?php echo $_GET['id']; ?>/mixed_train_balanced.csv','popup','width=1150,height=800'); return false;"><i class="fas fa-search-plus"></i> View</a>&emsp;<a href="/e/files/<?php echo $_GET['id']; ?>/mixed_train_balanced.csv" class="btn btn-warning"><i class="fas fa-download"></i> Download</a></td>
                        </tr>

                        <tr>
                            <td>Testing set:</td>
                            <td><a href="viewer.php?f=<?php echo $_GET['id']; ?>/mixed_test.csv" class="btn btn-info" role="button" target="popup" onclick="window.open('viewer.php?f=<?php echo $_GET['id']; ?>/mixed_test.csv','popup','width=1150,height=800'); return false;"><i class="fas fa-search-plus"></i> View</a>&emsp;<a href="/e/files/<?php echo $_GET['id']; ?>/mixed_test.csv" class="btn btn-warning"><i class="fas fa-download"></i> Download</a></td>
                        </tr>

                        <tr>
                            <td>Validation set:</td>
                            <td><a href="viewer.php?f=<?php echo $_GET['id']; ?>/mixed_valid.csv" class="btn btn-info" role="button" target="popup" onclick="window.open('viewer.php?f=<?php echo $_GET['id']; ?>/mixed_valid.csv','popup','width=1150,height=800'); return false;"><i class="fas fa-search-plus"></i> View</a>&emsp;<a href="/e/files/<?php echo $_GET['id']; ?>/mixed_valid.csv" class="btn btn-warning"><i class="fas fa-download"></i> Download</a></td>
                        </tr>

                        <tr>
                            <td>Merged data (training, testing, validation and balanced training):</td>
                            <td><a href="viewer.php?f=<?php echo $_GET['id']; ?>/merged.csv" class="btn btn-info" role="button" target="popup" onclick="window.open('viewer.php?f=<?php echo $_GET['id']; ?>/merged.csv','popup','width=1150,height=800'); return false;"><i class="fas fa-search-plus"></i> View</a>&emsp;<a href="/e/files/<?php echo $_GET['id']; ?>/merged.csv" class="btn btn-warning"><i class="fas fa-download"></i> Download</a></td>
                        </tr>

                        <tr>
                            <td>Differential expression analysis (training set only):</td>
                            <td><a href="viewer.php?f=<?php echo $_GET['id']; ?>/DE_train.csv" class="btn btn-info" role="button" target="popup" onclick="window.open('viewer.php?f=<?php echo $_GET['id']; ?>/DE_train.csv','popup','width=1150,height=800'); return false;"><i class="fas fa-search-plus"></i> View</a>&emsp;<a href="/e/files/<?php echo $_GET['id']; ?>/DE_train.csv" class="btn btn-warning"><i class="fas fa-download"></i> Download</a></td>
                        </tr>

                        <tr>
                            <td>Differential expression analysis (whole dataset):</td>
                            <td><a href="viewer.php?f=<?php echo $_GET['id']; ?>/DE_mixed.csv" class="btn btn-info" role="button" target="popup" onclick="window.open('viewer.php?f=<?php echo $_GET['id']; ?>/DE_mixed.csv','popup','width=1150,height=800'); return false;"><i class="fas fa-search-plus"></i> View</a>&emsp;<a href="/e/files/<?php echo $_GET['id']; ?>/DE_mixed.csv" class="btn btn-warning"><i class="fas fa-download"></i> Download</a></td>
                        </tr>

                        <tr>
                            <td>Re-do the initial check and default exploratory analysis:</td>
                            <td><a href="/e/notebooks/<?php echo $_GET['id']; ?>/formalcheckcsv.R" class="btn btn-danger" role="button" target="popup" onclick="window.open('/e/notebooks/<?php echo $_GET['id']; ?>/formalcheckcsv.R','popup','width=1150,height=800'); return false;"><i class="fas fa-play"></i> Run</a></td>
                        </tr>

                        <tr>
                            <td>Create your own analysis:</td>
                            <td><a href="/e/notebooks/<?php echo $_GET['id']; ?>/own_analysis.R" class="btn btn-danger" role="button" target="popup" onclick="window.open('/e/notebooks/<?php echo $_GET['id']; ?>/own_analysis.R','popup','width=1150,height=800'); return false;"><i class="fas fa-play"></i> Run</a></td>
                        </tr>

                    </tbody>
                </table>
            </div>
            <?php if ($var_initcheck != "OK") {
                die("There is an error in your data. See the log above. Fix this error, go back and try again with a new file. In case of any questions write to konrad.stawiski@umed.lodz.pl");
            }  ?>
            <?php if (!file_exists($target_dir . "featureselection_formulas_all.csv")) { ?>
                <div class="panel panel-primary">
                    <div class="panel-heading"><i class="fas fa-microscope"></i>&emsp;&emsp;Feature selection</div>
                    <div class="panel-body">
                        <p>
                            <font size="1">Notes: <i>By the defult all the methods are selected, but you can turn some of them off.</i></font>
                        </p>
                        <table class="table">
                            <form action="process.php?type=new_fs" method="post">
                                <input type="hidden" id="analysisid" name="analysisid" value="<?php echo $_GET['id']; ?>">
                                <thead>
                                    <th>Select</th>
                                    <th>ID</th>
                                    <th>Description</th>
                                    </td>
                                </thead>
                                <tbody>
                                    <tr>
                                        <td><label class="switch"><input type="checkbox" name="method[]" value="1" checked><span class="slider round"></span></label></td>
                                        <td>No: 1<br /><code>all<code></td>
                                        <td>Get all features (all features staring with 'hsa').</td>
                                    </tr>
                                    <tr>
                                        <td><label class="switch"><input type="checkbox" name="method[]" value="2" checked><span class="slider round"></span></label></td>
                                        <td>No: 2<br /><code>sig, sigtop, sigtopBonf, sigtopHolm, topFC, sigSMOTE, sigtopSMOTE, sigtopBonfSMOTE, sigtopHolmSMOTE, topFCSMOTE</code></td>
                                        <td>Selects features significantly differently expressed between classes by performing unpaired t-test with and without correction for multiple testing. We get: <code>sig</code> - all significant (adjusted p-value less or equal to 0.05) miRNAs with comparison using unpaired t-test and after the Benjamini-Hochberg procedure (BH, false discovery rate); <code>sigtop</code> - <code>sig</code> but limited only to the your prefered number of features (most significant features sorted by p-value), <code>sigtopBonf</code> - uses Bonferroni instead of BH correction, <code>sigtopHolm</code> - uses Holm–Bonferroni instead of BH correction, <code>topFC</code> - selects prefered number of features based on decreasing absolute value of fold change in differential analysis.
                                            <br />All the methods are also checked on dataset balanced with SMOTE (<a href="https://arxiv.org/pdf/1106.1813.pdf" target="_blank">Synthetic Minority Oversampling TEchnique</a>) - those formulas which names are appended with <code>SMOTE</code>.</td>
                                    </tr>

                                    <tr>
                                        <td><label class="switch"><input type="checkbox" name="method[]" value="3" checked><span class="slider round"></span></label></td>
                                        <td>No: 3<br /><code>fcsig, fcsigSMOTE</code></td>
                                        <td>Features significant in DE analysis using unpaired t-test and which absolute log2FC is greater than 1. Thus, features significant and up- or down-regulated in the higher magnitudes. FC - fold-change, DE - differential analysis.</td>
                                    </tr>

                                    <tr>
                                        <td><label class="switch"><input type="checkbox" name="method[]" value="4" checked><span class="slider round"></span></label></td>
                                        <td>No: 4<br /><code>cfs, cfsSMOTE, cfs_sig, cfsSMOTE_sig</code></td>
                                        <td><a href="https://www.cs.waikato.ac.nz/~mhall/thesis.pdf" target="_blank">Correlation-based feature selection</a> (CFS) - a heuristic algorithm selecting features that are highly correlated with class (binary) and lowly correlated with one another. It explores a search space in best-first manner, until stopping criteria are met.</td>
                                    </tr>


                                    <tr>
                                        <td><label class="switch"><input type="checkbox" name="method[]" value="5" checked><span class="slider round"></span></label></td>
                                        <td>No: 5<br /><code>classloop</code></td>
                                        <td>Classifier loop - performs multiple classification procedures using various algorithms (with embedded feature ranking) and various performance metrices. Final feature selection is done by combining the results. Modeling methods used: support vector machines, linear discriminant a nalysis, random forest and nearest shrunken centroid. Features are selected based on the AUC ROC and assessed in k-fold cross-validation according to the <a href="https://www.rdocumentation.org/packages/Biocomb/versions/0.4/topics/classifier.loop" target="_blank">documentation</a>. As this requires time, we do not perform it on SMOTEd dataset at once.</td>
                                    </tr>

                                    <tr>
                                        <td><label class="switch"><input type="checkbox" name="method[]" value="6" checked><span class="slider round"></span></label></td>
                                        <td>No: 6<br /><code>classloopSMOTE</code></td>
                                        <td>Application of <code>classloop</code> on balanced dataset (with SMOTE).</td>
                                    </tr>

                                    <tr>
                                        <td><label class="switch"><input type="checkbox" name="method[]" value="7" checked><span class="slider round"></span></label></td>
                                        <td>No: 7<br /><code>classloop_sig</code></td>
                                        <td>Application of <code>classloop</code> but only on the features which are significant in DE.</td>
                                    </tr>

                                    <tr>
                                        <td><label class="switch"><input type="checkbox" name="method[]" value="8" checked><span class="slider round"></span></label></td>
                                        <td>No: 8<br /><code>classloopSMOTE_sig</code></td>
                                        <td>Application of <code>classloop</code> on balanced training set and only on the features which are significant in DE (after balancing).</td>
                                    </tr>

                                    <tr>
                                        <td><label class="switch"><input type="checkbox" name="method[]" value="9" checked><span class="slider round"></span></label></td>
                                        <td>No: 9<br /><code>fcfs</code></td>
                                        <td>An algorithm similar to CFS, though exploring search space in greedy forward search manner (adding one, most attractive, feature at the time, until such addition does not improve set’s overall quality). Based on <a href="http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.411.9868&rep=rep1&type=pdf" target="_blank">Wang et al. 2005</a> and documented <a href="https://www.rdocumentation.org/packages/Biocomb/versions/0.4/topics/select.forward.Corr" target="_blank">here</a>.</td>
                                    </tr>

                                    <tr>
                                        <td><label class="switch"><input type="checkbox" name="method[]" value="10" checked><span class="slider round"></span></label></td>
                                        <td>No: 10<br /><code>fcfsSMOTE</code></td>
                                        <td>Application of <code>fcfs</code> on balanced training set.</td>
                                    </tr>

                                    <tr>
                                        <td><label class="switch"><input type="checkbox" name="method[]" value="11" checked><span class="slider round"></span></label></td>
                                        <td>No: 11<br /><code>fcfs_sig</code></td>
                                        <td>Application of <code>fcfs</code> on features significant in DE.</td>
                                    </tr>

                                    <tr>
                                        <td><label class="switch"><input type="checkbox" name="method[]" value="12" checked><span class="slider round"></span></label></td>
                                        <td>No: 12<br /><code>fcfsSMOTE_sig</code></td>
                                        <td>Application of <code>fcfs</code> on balanced dataset and on features significant in DE (after balancing).</td>
                                    </tr>

                                    <tr>
                                        <td><label class="switch"><input type="checkbox" name="method[]" value="13" checked><span class="slider round"></span></label></td>
                                        <td>No: 13<br /><code>fwrap</code></td>
                                        <td>A decision tree algorithm and forward search strategy documented <a href="https://www.rdocumentation.org/packages/Biocomb/versions/0.4/topics/select.forward.wrapper" target="_blank">here</a>.</td>
                                    </tr>

                                    <tr>
                                        <td><label class="switch"><input type="checkbox" name="method[]" value="14" checked><span class="slider round"></span></label></td>
                                        <td>No: 14<br /><code>fwrapSMOTE</code></td>
                                        <td>Application of <code>fwrap</code> on balanced training set.</td>
                                    </tr>

                                    <tr>
                                        <td><label class="switch"><input type="checkbox" name="method[]" value="15" checked><span class="slider round"></span></label></td>
                                        <td>No: 15<br /><code>fwrap_sig</code></td>
                                        <td>Application of <code>fwrap</code> on features significant in DE.</td>
                                    </tr>

                                    <tr>
                                        <td><label class="switch"><input type="checkbox" name="method[]" value="16" checked><span class="slider round"></span></label></td>
                                        <td>No: 16<br /><code>fwrapSMOTE_sig</code></td>
                                        <td>Application of <code>fwrap</code> on balanced dataset and on features significant in DE (after balancing).</td>
                                    </tr>

                                    <tr>
                                        <td><label class="switch"><input type="checkbox" name="method[]" value="17" checked><span class="slider round"></span></label></td>
                                        <td>No: 17<br /><code>AUC_MDL</code></td>
                                        <td>Feature ranking based on ROC AUC and minimal description length (MDL) discretization algorithm documented <a href="https://www.rdocumentation.org/packages/Biocomb/versions/0.4/topics/select.process" target="_blank">here</a>. After the ranking, the number of features are limited as set in options below.</td>
                                    </tr>

                                    <tr>
                                        <td><label class="switch"><input type="checkbox" name="method[]" value="18" checked><span class="slider round"></span></label></td>
                                        <td>No: 18<br /><code>SU_MDL</code></td>
                                        <td>Feature ranking based on symmetrical uncertainty and minimal description length (MDL) discretization algorithm documented <a href="https://www.rdocumentation.org/packages/Biocomb/versions/0.4/topics/select.process" target="_blank">here</a>. After the ranking, the number of features are limited as set in options below.</td>
                                    </tr>

                                    <tr>
                                        <td><label class="switch"><input type="checkbox" name="method[]" value="19" checked><span class="slider round"></span></label></td>
                                        <td>No: 19<br /><code>CorrSF_MDL</code></td>
                                        <td>Feature ranking based on CFS algorithm with forward search and minimal description length (MDL) discretization algorithm documented <a href="https://www.rdocumentation.org/packages/Biocomb/versions/0.4/topics/select.process" target="_blank">here</a>. After the ranking, the number of features are limited as set in options below.</td>
                                    </tr>

                                    <tr>
                                        <td><label class="switch"><input type="checkbox" name="method[]" value="20" checked><span class="slider round"></span></label></td>
                                        <td>No: 20<br /><code>AUC_MDLSMOTE</code></td>
                                        <td>Feature ranking based on ROC AUC and minimal description length (MDL) discretization algorithm documented <a href="https://www.rdocumentation.org/packages/Biocomb/versions/0.4/topics/select.process" target="_blank">here</a>. After the ranking, the number of features are limited as set in options below. Performed on the training set balanced with SMOTE.</td>
                                    </tr>

                                    <tr>
                                        <td><label class="switch"><input type="checkbox" name="method[]" value="21" checked><span class="slider round"></span></label></td>
                                        <td>No: 21<br /><code>SU_MDLSMOTE</code></td>
                                        <td>Feature ranking based on symmetrical uncertainty and minimal description length (MDL) discretization algorithm documented <a href="https://www.rdocumentation.org/packages/Biocomb/versions/0.4/topics/select.process" target="_blank">here</a>. After the ranking, the number of features are limited as set in options below. Performed on the training set balanced with SMOTE.</td>
                                    </tr>

                                    <tr>
                                        <td><label class="switch"><input type="checkbox" name="method[]" value="22" checked><span class="slider round"></span></label></td>
                                        <td>No: 22<br /><code>CorrSF_MDLSMOTE</code></td>
                                        <td>Feature ranking based on CFS algorithm with forward search and minimal description length (MDL) discretization algorithm documented <a href="https://www.rdocumentation.org/packages/Biocomb/versions/0.4/topics/select.process" target="_blank">here</a>. After the ranking, the number of features are limited as set in options below. Performed on the training set balanced with SMOTE.</td>
                                    </tr>

                                    <tr>
                                        <td><label class="switch"><input type="checkbox" name="method[]" value="23" checked><span class="slider round"></span></label></td>
                                        <td>No: 23<br /><code>AUC_MDL_sig</code></td>
                                        <td>Feature ranking based on ROC AUC and minimal description length (MDL) discretization algorithm documented <a href="https://www.rdocumentation.org/packages/Biocomb/versions/0.4/topics/select.process" target="_blank">here</a>. After the ranking, the number of features are limited as set in options below. Only features significant in DE are allowed.</td>
                                    </tr>

                                    <tr>
                                        <td><label class="switch"><input type="checkbox" name="method[]" value="24" checked><span class="slider round"></span></label></td>
                                        <td>No: 24<br /><code>SU_MDL_sig</code></td>
                                        <td>Feature ranking based on symmetrical uncertainty and minimal description length (MDL) discretization algorithm documented <a href="https://www.rdocumentation.org/packages/Biocomb/versions/0.4/topics/select.process" target="_blank">here</a>. After the ranking, the number of features are limited as set in options below. Only features significant in DE are allowed.</td>
                                    </tr>

                                    <tr>
                                        <td><label class="switch"><input type="checkbox" name="method[]" value="25" checked><span class="slider round"></span></label></td>
                                        <td>No: 25<br /><code>CorrSF_MDL_sig</code></td>
                                        <td>Feature ranking based on CFS algorithm with forward search and minimal description length (MDL) discretization algorithm documented <a href="https://www.rdocumentation.org/packages/Biocomb/versions/0.4/topics/select.process" target="_blank">here</a>. After the ranking, the number of features are limited as set in options below. Only features significant in DE are allowed.</td>
                                    </tr>

                                    <tr>
                                        <td><label class="switch"><input type="checkbox" name="method[]" value="26" checked><span class="slider round"></span></label></td>
                                        <td>No: 26<br /><code>AUC_MDLSMOTE_sig</code></td>
                                        <td>Feature ranking based on ROC AUC and minimal description length (MDL) discretization algorithm documented <a href="https://www.rdocumentation.org/packages/Biocomb/versions/0.4/topics/select.process" target="_blank">here</a>. After the ranking, the number of features are limited as set in options below. Performed on the training set balanced with SMOTE. Only features significant in DE are allowed.</td>
                                    </tr>

                                    <tr>
                                        <td><label class="switch"><input type="checkbox" name="method[]" value="27" checked><span class="slider round"></span></label></td>
                                        <td>No: 27<br /><code>SU_MDLSMOTE_sig</code></td>
                                        <td>Feature ranking based on symmetrical uncertainty and minimal description length (MDL) discretization algorithm documented <a href="https://www.rdocumentation.org/packages/Biocomb/versions/0.4/topics/select.process" target="_blank">here</a>. After the ranking, the number of features are limited as set in options below. Performed on the training set balanced with SMOTE. Only features significant in DE are allowed.</td>
                                    </tr>

                                    <tr>
                                        <td><label class="switch"><input type="checkbox" name="method[]" value="28" checked><span class="slider round"></span></label></td>
                                        <td>No: 28<br /><code>CorrSF_MDLSMOTE_sig</code></td>
                                        <td>Feature ranking based on CFS algorithm with forward search and minimal description length (MDL) discretization algorithm documented <a href="https://www.rdocumentation.org/packages/Biocomb/versions/0.4/topics/select.process" target="_blank">here</a>. After the ranking, the number of features are limited as set in options below. Performed on the training set balanced with SMOTE. Only features significant in DE are allowed.</td>
                                    </tr>

                                    <tr>
                                        <td><label class="switch"><input type="checkbox" name="method[]" value="29" checked><span class="slider round"></span></label></td>
                                        <td>No: 29<br /><code>bounceR-full, bounceR-stability</code></td>
                                        <td>A component-wise-boosting-based algorithm selecting optimal features in multiple iterations of single feature-models construction. See the source <a href="https://github.com/STATWORX/bounceR" target="_blank">here</a>. <code>bounceR-stability</code> gets the most stable features. Wrapper methods implemented here leverage componentwise boosting as a weak learners.</td>
                                    </tr>

                                    <tr>
                                        <td><label class="switch"><input type="checkbox" name="method[]" value="30"><span class="slider round"></span></label></td>
                                        <td>No: 30<br /><code>bounceR-full_SMOTE, bounceR-stability_SMOTE</code></td>
                                        <td>A component-wise-boosting-based algorithm selecting optimal features in multiple iterations of single feature-models construction. See the source <a href="https://github.com/STATWORX/bounceR" target="_blank">here</a>. <code>bounceR-stability</code> gets the most stable features. Wrapper methods implemented here leverage componentwise boosting as a weak learners. Performed on the training set balanced with SMOTE. <br /><i>Note: This method may take a lot of memory or other resources for larger datasets, thus is disabled by default!</i></td>
                                    </tr>

                                    <tr>
                                        <td><label class="switch"><input type="checkbox" name="method[]" value="31" checked><span class="slider round"></span></label></td>
                                        <td>No: 31<br /><code>bounceR-full_SIG, bounceR-stability_SIG</code></td>
                                        <td>A component-wise-boosting-based algorithm selecting optimal features in multiple iterations of single feature-models construction. See the source <a href="https://github.com/STATWORX/bounceR" target="_blank">here</a>. <code>bounceR-stability</code> gets the most stable features. Wrapper methods implemented here leverage componentwise boosting as a weak learners. Only features significant in DE are allowed. </td>
                                    </tr>

                                    <tr>
                                        <td><label class="switch"><input type="checkbox" name="method[]" value="32"><span class="slider round"></span></label></td>
                                        <td>No: 32<br /><code>bounceR-full_SIGSMOTE, bounceR-stability_SIGSMOTE</code></td>
                                        <td>A component-wise-boosting-based algorithm selecting optimal features in multiple iterations of single feature-models construction. See the source <a href="https://github.com/STATWORX/bounceR" target="_blank">here</a>. <code>bounceR-stability</code> gets the most stable features. Wrapper methods implemented here leverage componentwise boosting as a weak learners. Only features significant in DE are allowed. Performed on the training set balanced with SMOTE. <br /><i>Note: This method may take a lot of memory or other resources for larger datasets, thus is disabled by default!</i> </td>
                                    </tr>

                                    <tr>
                                        <td><label class="switch"><input type="checkbox" name="method[]" value="33" checked><span class="slider round"></span></label></td>
                                        <td>No: 33<br /><code>RandomForestRFE</code></td>
                                        <td>Recursively eliminates features from the feature space based on ranking from Random Forrest classifier (retrained woth resampling after each elimination). Details are available <a href="https://topepo.github.io/caret/recursive-feature-elimination.html#search" target="_blank">here</a>. </td>
                                    </tr>

                                    <tr>
                                        <td><label class="switch"><input type="checkbox" name="method[]" value="34"><span class="slider round"></span></label></td>
                                        <td>No: 34<br /><code>RandomForestRFESMOTE</code></td>
                                        <td>Recursively eliminates features from the feature space based on ranking from Random Forrest classifier (retrained woth resampling after each elimination). Details are available <a href="https://topepo.github.io/caret/recursive-feature-elimination.html#search" target="_blank">here</a>. Performed on the training set balanced with SMOTE.<br /><i>Note: This method may take a lot of memory or other resources for larger datasets, thus is disabled by default!</i> </td>
                                    </tr>

                                    <tr>
                                        <td><label class="switch"><input type="checkbox" name="method[]" value="35" checked><span class="slider round"></span></label></td>
                                        <td>No: 35<br /><code>RandomForestRFE_sig</code></td>
                                        <td>Recursively eliminates features from the feature space based on ranking from Random Forrest classifier (retrained woth resampling after each elimination). Details are available <a href="https://topepo.github.io/caret/recursive-feature-elimination.html#search" target="_blank">here</a>. Only features significant in DE are allowed. </td>
                                    </tr>

                                    <tr>
                                        <td><label class="switch"><input type="checkbox" name="method[]" value="36"><span class="slider round"></span></label></td>
                                        <td>No: 36<br /><code>RandomForestRFESMOTE_sig</code></td>
                                        <td>Recursively eliminates features from the feature space based on ranking from Random Forrest classifier (retrained woth resampling after each elimination). Details are available <a href="https://topepo.github.io/caret/recursive-feature-elimination.html#search" target="_blank">here</a>. Only features significant in DE are allowed. Performed on the training set balanced with SMOTE. <br /><i>Note: This method may take a lot of memory or other resources for larger datasets, thus is disabled by default!</i> </td>
                                    </tr>

                                    <tr>
                                        <td><label class="switch"><input type="checkbox" name="method[]" value="37" checked><span class="slider round"></span></label></td>
                                        <td>No: 37<br /><code>GeneticAlgorithmRF</code></td>
                                        <td>Uses genetic algorithm principle to search for optimal subset of the feature space. This uses internally implemented random forest model and 10-fold cross validation to assess performance of the "chromosomes" in each generation. Details are available <a href="https://topepo.github.io/caret/feature-selection-using-genetic-algorithms.html" target="_blank">here</a>. </td>
                                    </tr>

                                    <tr>
                                        <td><label class="switch"><input type="checkbox" name="method[]" value="38"><span class="slider round"></span></label></td>
                                        <td>No: 38<br /><code>GeneticAlgorithmRFSMOTE</code></td>
                                        <td>Uses genetic algorithm principle to search for optimal subset of the feature space. This uses internally implemented random forest model and 10-fold cross validation to assess performance of the "chromosomes" in each generation. Details are available <a href="https://topepo.github.io/caret/feature-selection-using-genetic-algorithms.html" target="_blank">here</a>. Performed on the training set balanced with SMOTE. <br /><i>Note: This method may take a lot of memory or other resources for larger datasets, thus is disabled by default!</i> </td>
                                    </tr>

                                    <tr>
                                        <td><label class="switch"><input type="checkbox" name="method[]" value="39" checked><span class="slider round"></span></label></td>
                                        <td>No: 39<br /><code>GeneticAlgorithmRF_sig</code></td>
                                        <td>Uses genetic algorithm principle to search for optimal subset of the feature space. This uses internally implemented random forest model and 10-fold cross validation to assess performance of the "chromosomes" in each generation. Details are available <a href="https://topepo.github.io/caret/feature-selection-using-genetic-algorithms.html" target="_blank">here</a>. Only features significant in DE are allowed. </td>
                                    </tr>

                                    <tr>
                                        <td><label class="switch"><input type="checkbox" name="method[]" value="40"><span class="slider round"></span></label></td>
                                        <td>No: 40<br /><code>GeneticAlgorithmRFSMOTE_sig</code></td>
                                        <td>Uses genetic algorithm principle to search for optimal subset of the feature space. This uses internally implemented random forest model and 10-fold cross validation to assess performance of the "chromosomes" in each generation. Details are available <a href="https://topepo.github.io/caret/feature-selection-using-genetic-algorithms.html" target="_blank">here</a>. Only features significant in DE are allowed. Performed on the training set balanced with SMOTE. <br /><i>Note: This method may take a lot of memory or other resources for larger datasets, thus is disabled by default!</i> </td>
                                    </tr>

                                    <tr>
                                        <td><label class="switch"><input type="checkbox" name="method[]" value="41" checked><span class="slider round"></span></label></td>
                                        <td>No: 41<br /><code>SimulatedAnnealingRF</code></td>
                                        <td>Simulated Annealing - explores a feature space by randomly modifying a given feature subset and evaluating classification performance using new attributes to check whether changes were beneficial. It is is a global search method that makes small random changes (i.e. perturbations) to an initial candidate solution. In this method also random forest is used as a model for evaluation. Details are available <a href="https://topepo.github.io/caret/feature-selection-using-simulated-annealing.html" target="_blank">here</a>. </td>
                                    </tr>

                                    <tr>
                                        <td><label class="switch"><input type="checkbox" name="method[]" value="42" checked><span class="slider round"></span></label></td>
                                        <td>No: 42<br /><code>SimulatedAnnealingRFSMOTE</code></td>
                                        <td>Simulated Annealing - explores a feature space by randomly modifying a given feature subset and evaluating classification performance using new attributes to check whether changes were beneficial. It is is a global search method that makes small random changes (i.e. perturbations) to an initial candidate solution. In this method also random forest is used as a model for evaluation. Details are available <a href="https://topepo.github.io/caret/feature-selection-using-simulated-annealing.html" target="_blank">here</a>. Performed on the training set balanced with SMOTE. </td>
                                    </tr>

                                    <tr>
                                        <td><label class="switch"><input type="checkbox" name="method[]" value="43" checked><span class="slider round"></span></label></td>
                                        <td>No: 43<br /><code>SimulatedAnnealingRF_sig</code></td>
                                        <td>Simulated Annealing - explores a feature space by randomly modifying a given feature subset and evaluating classification performance using new attributes to check whether changes were beneficial. It is is a global search method that makes small random changes (i.e. perturbations) to an initial candidate solution. In this method also random forest is used as a model for evaluation. Details are available <a href="https://topepo.github.io/caret/feature-selection-using-simulated-annealing.html" target="_blank">here</a>. Only features significant in DE are allowed. </td>
                                    </tr>

                                    <tr>
                                        <td><label class="switch"><input type="checkbox" name="method[]" value="44" checked><span class="slider round"></span></label></td>
                                        <td>No: 44<br /><code>SimulatedAnnealingRFSMOTE_sig</code></td>
                                        <td>Simulated Annealing - explores a feature space by randomly modifying a given feature subset and evaluating classification performance using new attributes to check whether changes were beneficial. It is is a global search method that makes small random changes (i.e. perturbations) to an initial candidate solution. In this method also random forest is used as a model for evaluation. Details are available <a href="https://topepo.github.io/caret/feature-selection-using-simulated-annealing.html" target="_blank">here</a>. Only features significant in DE are allowed. Performed on the training set balanced with SMOTE. </td>
                                    </tr>

                                    <tr>
                                        <td><label class="switch"><input type="checkbox" name="method[]" value="45" checked><span class="slider round"></span></label></td>
                                        <td>No: 45<br /><code>Boruta</code></td>
                                        <td>Boruta - utilizes random forrest algorithm to iteratively remove features proved to be less relevant than random variables. Details are available in paper by <a href="https://www.jstatsoft.org/v36/i11/paper/" target="_blank">Kursa et al. 2010</a> or <a href="https://www.datacamp.com/community/tutorials/feature-selection-R-boruta" target="_blank">this blog post</a>. Only features significant in DE are allowed. Performed on the training set balanced with SMOTE. </td>
                                    </tr>

                                    <tr>
                                        <td><label class="switch"><input type="checkbox" name="method[]" value="46" checked><span class="slider round"></span></label></td>
                                        <td>No: 46<br /><code>BorutaSMOTE</code></td>
                                        <td>Boruta - utilizes random forrest algorithm to iteratively remove features proved to be less relevant than random variables. Details are available in paper by <a href="https://www.jstatsoft.org/v36/i11/paper/" target="_blank">Kursa et al. 2010</a> or <a href="https://www.datacamp.com/community/tutorials/feature-selection-R-boruta" target="_blank">this blog post</a>. Performed on the training set balanced with SMOTE. </td>
                                    </tr>

                                    <tr>
                                        <td><label class="switch"><input type="checkbox" name="method[]" value="47" checked><span class="slider round"></span></label></td>
                                        <td>No: 47<br /><code>spFSR</code></td>
                                        <td>spFSR - feature selection and ranking by simultaneous perturbation stochastic approximation. This is an algorithm based on pseudo-gradient descent stochastic optimisation with Barzilai-Borwein method for step size and gradient estimation optimization. Details are available in paper by <a href="https://arxiv.org/abs/1804.05589" target="_blank">Zeren et al. 2018</a>. </td>
                                    </tr>

                                    <tr>
                                        <td><label class="switch"><input type="checkbox" name="method[]" value="48" checked><span class="slider round"></span></label></td>
                                        <td>No: 48<br /><code>spFSRSMOTE</code></td>
                                        <td>spFSR - feature selection and ranking by simultaneous perturbation stochastic approximation. This is an algorithm based on pseudo-gradient descent stochastic optimisation with Barzilai-Borwein method for step size and gradient estimation optimization. Details are available in paper by <a href="https://arxiv.org/abs/1804.05589" target="_blank">Zeren et al. 2018</a>. Performed on the training set balanced with SMOTE. </td>
                                    </tr>

                                    <tr>
                                        <td><label class="switch"><input type="checkbox" name="method[]" value="49" checked><span class="slider round"></span></label></td>
                                        <td>No: 49<br /><code>varSelRF, varSelRFSMOTE</code></td>
                                        <td>varSelRF - recursively eliminates features using random forrest feature scores, seeking to minimize out-of-bag classification error. Details are available in paper by <a href="https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-7-3" target="_blank">Díaz-Uriarte et al. 2006</a>. Performed on the unbalanced training set as well as on balanced with SMOTE. </td>
                                    </tr>

                                    <tr>
                                        <td><label class="switch"><input type="checkbox" name="method[]" value="50" checked><span class="slider round"></span></label></td>
                                        <td>No: 50<br /><code>Wx, WxSMOTE</code></td>
                                        <td>Wx - deep neural network-based (deep learning) feature (gene) selection algorithm. <a href="https://github.com/kstawiski/OmicSelector/blob/master/inst/extdata/wx/DearWXpub/src/wx_konsta.py" target="_blank">We use 2 hidden layers with 16 hidden neurons.</a> Details are available in paper by <a href="https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6642261/" target="_blank">Park et al. 2019</a>. Performed on the unbalanced training set as well as on balanced with SMOTE. </td>
                                    </tr>

                                    <tr>
                                        <td><label class="switch"><input type="checkbox" name="method[]" value="51" checked><span class="slider round"></span></label></td>
                                        <td>No: 51<br /><code>Mystepwise_glm_binomial, Mystepwise_sig_glm_binomial</code></td>
                                        <td>Stepwise variable selection procedure (with iterations between the 'forward' and 'backward' steps) for generalized linear models with logit link function (i.e. logistic regression). We use p=0.05 as a threshold for both entry (SLE) and stay (SLS). Details are available <a href="https://www.rdocumentation.org/packages/My.stepwise" target="_blank">here</a>. Performed on all features of training set as well as features initially selected in DE (significant in DE). </td>
                                    </tr>

                                    <tr>
                                        <td><label class="switch"><input type="checkbox" name="method[]" value="51" checked><span class="slider round"></span></label></td>
                                        <td>No: 52<br /><code>Mystepwise_glm_binomialSMOTE, Mystepwise_sig_glm_binomialSMOTE</code></td>
                                        <td>Stepwise variable selection procedure (with iterations between the 'forward' and 'backward' steps) for generalized linear models with logit link function (i.e. logistic regression). We use p=0.05 as a threshold for both entry (SLE) and stay (SLS). Details are available <a href="https://www.rdocumentation.org/packages/My.stepwise" target="_blank">here</a>. Performed on all features of training set as well as features initially selected in DE (significant in DE) after balancing the training set with SMOTE. </td>
                                    </tr>

                                    <tr>
                                        <td><label class="switch"><input type="checkbox" name="method[]" value="53" checked><span class="slider round"></span></label></td>
                                        <td>No: 53<br /><code>stepAIC, stepAICsig</code></td>
                                        <td>Here we perform a stepwise model selection by AIC (Akaike Information Criterion) based on logistic regression. Details are available <a href="https://www.rdocumentation.org/packages/MASS/versions/7.3-53/topics/stepAIC" target="_blank">here</a>. Performed on all features of training set as well as features initially selected in DE (significant in DE). </td>
                                    </tr>

                                    <tr>
                                        <td><label class="switch"><input type="checkbox" name="method[]" value="54" checked><span class="slider round"></span></label></td>
                                        <td>No: 54<br /><code>stepAIC_SMOTE, stepAICsig_SMOTE</code></td>
                                        <td>Here we perform a stepwise model selection by AIC (Akaike Information Criterion) based on logistic regression. Details are available <a href="https://www.rdocumentation.org/packages/MASS/versions/7.3-53/topics/stepAIC" target="_blank">here</a>. Performed on all features of training set as well as features initially selected in DE (significant in DE) after balancing the training set with SMOTE. </td>
                                    </tr>

                                    <tr>
                                        <td><label class="switch"><input type="checkbox" name="method[]" value="55"><span class="slider round"></span></label></td>
                                        <td>No: 55<br /><code>iteratedRFECV, iteratedRFETest</code></td>
                                        <td>Iterated RFE tested in cross-validation and on test set (watch out for bias!). See the source <a href="https://github.com/kstawiski/OmicSelector/blob/master/R/OmicSelector_iteratedRFE.R" target="_blank">here</a>. <br /><i>Note: This method may take a lot of memory or other resources for larger datasets, thus is disabled by default!</i></td>
                                    </tr>

                                    <tr>
                                        <td><label class="switch"><input type="checkbox" name="method[]" value="56"><span class="slider round"></span></label></td>
                                        <td>No: 56<br /><code>iteratedRFECV_SMOTE, iteratedRFETest_SMOTE</code></td>
                                        <td>Iterated RFE tested in cross-validation and on test set (watch out for bias!). See the source <a href="https://github.com/kstawiski/OmicSelector/blob/master/R/OmicSelector_iteratedRFE.R" target="_blank">here</a>. Performed after balancing the training set with SMOTE. <br /><i>Note: This method may take a lot of memory or other resources for larger datasets, thus is disabled by default!</i> </td>
                                    </tr>

                                    <tr>
                                        <td><label class="switch"><input type="checkbox" name="method[]" value="57" checked><span class="slider round"></span></label></td>
                                        <td>No: 57<br /><code>LASSO, LASSO_SMOTE</code></td>
                                        <td>Feature selection based on LASSO (Least Absolute Shrinkage and Selection Operator) model with alpha = 1 - penalizes with L1-norm; with 10-fold cross-validation. See the source <a href="https://www.rdocumentation.org/packages/glmnet/versions/4.0-2/topics/glmnet" target="_blank">here</a>. Performed on originial training set and after balancing the training set with SMOTE. </td>
                                    </tr>

                                    <tr>
                                        <td><label class="switch"><input type="checkbox" name="method[]" value="58" checked><span class="slider round"></span></label></td>
                                        <td>No: 58<br /><code>ElasticNet, ElasticNet_SMOTE</code></td>
                                        <td>Feature selection based on elastic net with tuning the value of alpha through a line search. See the source <a href="https://www.rdocumentation.org/packages/glmnet/versions/4.0-2/topics/glmnet" target="_blank">here</a>. Performed on originial training set and after balancing the training set with SMOTE. </td>
                                    </tr>

                                    <tr>
                                        <td><label class="switch"><input type="checkbox" name="method[]" value="59" checked><span class="slider round"></span></label></td>
                                        <td>No: 59<br /><code>stepLDA, stepLDA_SMOTE</code></td>
                                        <td>Forward/backward variable selection (both directions) for linear discriminant analysis. See the source <a href="https://www.rdocumentation.org/packages/klaR" target="_blank">here</a>. Performed on originial training set and after balancing the training set with SMOTE. </td>
                                    </tr>

                                    <tr>
                                        <td><label class="switch"><input type="checkbox" name="method[]" value="60" checked><span class="slider round"></span></label></td>
                                        <td>No: 60<br /><code>feseR_filter.corr, feseR_gain.inf, feseR_matrix.corr, feseR_combineFS_RF, feseR_filter.corr_SMOTE, feseR_gain.inf_SMOTE, feseR_matrix.corr_SMOTE, feseR_combineFS_RF_SMOTE</code></td>
                                        <td>Set of feature selection methods embeded in feseR package published by Perez-Rivelor et al. All default parameters are used, but mincorr is set to 0.2. See the paper <a href="https://doi.org/10.1371/journal.pone.0189875" target="_blank">here</a>. Performed on originial training set and after balancing the training set with SMOTE. </td>
                                    </tr>

                                    <tr>
                                        <td>Options:</td>
                                        <td colspan="2">
                                            <div class="form-group row">
                                                <div class="col-sm-7">
                                                    <p><u>Prefered number of features</u><br />
                                                        <font size="1"><i>(some of the methods do not select features but rank them, how many features are acceptable for you?)</i></font>
                                                    </p>
                                                </div>
                                                <div class="col-sm-5">
                                                    <input class="form-control" id="prefer_no_features" name="prefer_no_features" type="number" min="1" max="50" value="10" />
                                                </div>
                                            </div>
                                            <div class="form-group row">
                                                <div class="col-sm-7">
                                                    <p><u>Timeout for selected methods</u><br />
                                                        <font size="1"><i>(max time for the method to run in seconds, if you do not want to wait the ethernity for the results in misconfigured pipeline)</i></font>
                                                    </p>
                                                </div>
                                                <div class="col-sm-5">
                                                    <input class="form-control" id="timeout_sec" name="timeout_sec" type="number" min="0" max="2629743" value="3600" />
                                                </div>
                                            </div>
                                            <div class="form-group row">
                                                <div class="col-sm-7">
                                                    <p><u>Maximum number of iterations</u><br />
                                                        <font size="1"><i>(maximum number of iterations in selected methods, setting this too high may results in very long comupting time)</i></font>
                                                    </p>
                                                </div>
                                                <div class="col-sm-5">
                                                    <input class="form-control" id="max_iterations" name="max_iterations" type="number" min="1" max="50" value="5" />
                                                </div>
                                            </div>
                                        </td>
                                    </tr>
                                </tbody>
                        </table>
                        <p>
                            <button type="submit" class="btn btn-success" value="Upload" name="submit" onclick="waitingDialog.show('Starting the analysis...');">
                                <i class="fas fa-clipboard-check"></i>&emsp;Start feature selection
                            </button>&emsp;<a href="process.php?type=recover_fs&analysisid=<?php echo $_GET['id']; ?>" class="btn btn-primary" role="button"><i class="fas fa-trash-restore"></i> Recover interrupted selection</a></p>
                        </form>
                    </div>
                </div>
            <?php } else { ?>
                <div class="panel panel-success">
                    <div class="panel-heading"><i class="fas fa-microscope"></i>&emsp;&emsp;Feature selection</div>
                    <div class="panel-body">
                        <p>
                            <h3>Final set of feature sets selected for further evaluation:</h3>
                            <br />
                            <font size="1">Notes: This table presents final formulas. <a href="https://biostat.umed.pl/OmicSelector/reference/OmicSelector_merge_formulas.html" target="_blank">The formulas with features more than the prefered number of features of features were trimmed (according to documentation).</a> Using the download button you can download filtered set.</font><br>
                        </p>
                        <p><?php konsta_readcsv_formulas($target_dir . "featureselection_formulas_final.csv"); ?></p>
                        <p>
                            <h4>Details:</h4>

                            <table class="table">
                                <tr>
                                    <td>All formulas selected by the methods:</td>
                                    <td><a href="viewer.php?f=<?php echo $_GET['id']; ?>/featureselection_formulas_all.csv" class="btn btn-info" role="button" target="popup" onclick="window.open('viewer.php?f=<?php echo $_GET['id']; ?>/featureselection_formulas_all.csv','popup','width=1150,height=800'); return false;"><i class="fas fa-search-plus"></i> View</a>&emsp;<a href="/e/files/<?php echo $_GET['id']; ?>/featureselection_formulas_all.csv" class="btn btn-warning"><i class="fas fa-download"></i> Download</a></td>
                                </tr>
                                <tr>
                                    <td>Final set of formulas selected by the methods:<br />
                                        <font size="1"><i>(formulas with more than prefered number of features were removed, classical selection methods are intact)</i></font>
                                    <td><a href="viewer.php?f=<?php echo $_GET['id']; ?>/featureselection_formulas_final.csv" class="btn btn-info" role="button" target="popup" onclick="window.open('viewer.php?f=<?php echo $_GET['id']; ?>/featureselection_formulas_final.csv','popup','width=1150,height=800'); return false;"><i class="fas fa-search-plus"></i> View</a>&emsp;<a href="/e/files/<?php echo $_GET['id']; ?>/featureselection_formulas_final.csv" class="btn btn-warning"><i class="fas fa-download"></i> Download</a></td>
                                </tr>
                                <tr>
                                    <td>Options:</td>
                                    <td>
                                        <p><a href="/e/edit/<?php echo $_GET['id']; ?>/temp/featureselection.log" class="btn btn-primary" role="button" target="popup" onclick="window.open('/e/edit/<?php echo $_GET['id']; ?>/temp/featureselection.log','popup','width=1150,height=800'); return false;"><i class="fas fa-history"></i> View log</a>&emsp;<a href="process.php?type=delete_fs2&analysisid=<?php echo $_GET['id']; ?>" class="btn btn-success" role="button"><i class="fas fa-plus"></i> Add or re-run some methods</a></p>
                                        <p><a href="process.php?type=delete_fs&analysisid=<?php echo $_GET['id']; ?>" class="btn btn-danger" role="button" onclick="return confirm('Are you sure? This will delete all the data regarding your feature selection')"><i class="fas fa-trash"></i> Re-run selection (delete previous)</a></p>


                                    </td>
                                </tr>

                            </table>
                        </p>

                    </div>
                </div>
                <?php
                // Po selekcji zmiennych, przygotowanie benchmarku.
                if (!file_exists($target_dir . "benchmark.csv")) { ?>
                    <div class="panel panel-primary">
                        <div class="panel-heading"><i class="fas fa-award"></i>&emsp;&emsp;Benchmarking</div>
                        <div class="panel-body">
                            <p>
                                <font size="1">Notes: <i>By the defult all the methods are selected, but you can turn some of them off.</i></font>
                            </p>
                            <table class="table">
                                <form action="process.php?type=new_benchmark" method="post">
                                    <input type="hidden" id="analysisid" name="analysisid" value="<?php echo $_GET['id']; ?>">
                                    <thead>
                                        <tr>
                                            <th>Select</th>
                                            <th>ID</th>
                                            <th>Description</th>
                                        </tr>
                                    </thead>
                                    <tbody>
                                        <tr>
                                            <td>
                                                <font size="1">Always enabled.</font>
                                            </td>
                                            <td><code>glm<code></td>
                                            <td>Logistic regression (generalized linear model with binomial link function).</td>
                                        </tr>
                                        <tr>
                                            <td><label class="switch"><input type="checkbox" name="method[]" value="mlp" checked><span class="slider round"></span></label></td>
                                            <td><code>mlp<code></td>
                                            <td>Multilayer perceptron (MLP) - fully connected feedforward neural network with 1 hidden layer and logistic activiation function. Details: <a href="https://github.com/topepo/caret/blob/master/models/files/mlp.R" target="_blank">code</a>, <a href="https://www.rdocumentation.org/packages/RSNNS" target="_blank">package</a>.</td>
                                        </tr>
                                        <tr>
                                            <td><label class="switch"><input type="checkbox" name="method[]" value="mlpML" checked><span class="slider round"></span></label></td>
                                            <td><code>mlpML<code></td>
                                            <td>Multilayer perceptron (MLP) - fully connected feedforward neural network with up to 3 hidden layers and logistic activiation function. Details: <a href="https://github.com/topepo/caret/blob/master/models/files/mlpML.R" target="_blank">code</a>, <a href="https://www.rdocumentation.org/packages/RSNNS" target="_blank">package</a>.</td>
                                        </tr>
                                        <tr>
                                            <td><label class="switch"><input type="checkbox" name="method[]" value="svmRadial" checked><span class="slider round"></span></label></td>
                                            <td><code>svmRadial<code></td>
                                            <td>Support vector machines with radial basis function kernel. Details: <a href="https://github.com/topepo/caret/blob/master/models/files/svmRadial.R" target="_blank">code</a>, <a href="https://www.rdocumentation.org/packages/kernlab" target="_blank">package</a>.</td>
                                        </tr>
                                        <tr>
                                            <td><label class="switch"><input type="checkbox" name="method[]" value="svmLinear" checked><span class="slider round"></span></label></td>
                                            <td><code>svmLinear<code></td>
                                            <td>Support vector machines with linear kernel. Details: <a href="https://github.com/topepo/caret/blob/master/models/files/svmLinear.R" target="_blank">code</a>, <a href="https://www.rdocumentation.org/packages/kernlab" target="_blank">package</a>.</td>
                                        </tr>
                                        <tr>
                                            <td><label class="switch"><input type="checkbox" name="method[]" value="rf" checked><span class="slider round"></span></label></td>
                                            <td><code>rf<code></td>
                                            <td>Random forest. Details: <a href="https://github.com/topepo/caret/blob/master/models/files/rf.R" target="_blank">code</a>, <a href="https://www.rdocumentation.org/packages/randomForest" target="_blank">package</a>.</td>
                                        </tr>
                                        <tr>
                                            <td><label class="switch"><input type="checkbox" name="method[]" value="C5.0" checked><span class="slider round"></span></label></td>
                                            <td><code>C5.0<code></td>
                                            <td>C5.0 decision trees and rule-based models. Details: <a href="https://github.com/topepo/caret/blob/master/models/files/C5.0.R" target="_blank">code</a>, <a href="https://www.rdocumentation.org/packages/C50" target="_blank">package</a>.</td>
                                        </tr>
                                        <tr>
                                            <td><label class="switch"><input type="checkbox" name="method[]" value="rpart" checked><span class="slider round"></span></label></td>
                                            <td><code>rpart<code></td>
                                            <td>CART decision trees with modulation of complexity parameter. Details: <a href="https://github.com/topepo/caret/blob/master/models/files/rpart.R" target="_blank">code</a>, <a href="https://www.rdocumentation.org/packages/rpart" target="_blank">package</a>.</td>
                                        </tr>
                                        <tr>
                                            <td><label class="switch"><input type="checkbox" name="method[]" value="rpart2" checked><span class="slider round"></span></label></td>
                                            <td><code>rpart2<code></td>
                                            <td>CART decision trees with modulation of max tree depth. Details: <a href="https://github.com/topepo/caret/blob/master/models/files/rpart2.R" target="_blank">code</a>, <a href="https://www.rdocumentation.org/packages/rpart" target="_blank">package</a>.</td>
                                        </tr>
                                        <tr>
                                            <td><label class="switch"><input type="checkbox" name="method[]" value="ctree" checked><span class="slider round"></span></label></td>
                                            <td><code>ctree<code></td>
                                            <td>Conditional inference trees. Details: <a href="https://github.com/topepo/caret/blob/master/models/files/ctree.R" target="_blank">code</a>, <a href="https://www.rdocumentation.org/packages/party" target="_blank">package</a>.</td>
                                        </tr>
                                        <!-- <tr>
                                            <td><label class="switch"><input type="checkbox" name="mxnet" value="TRUE"><span class="slider round"></span></label></td>
                                            <td><code>mxnet<code></td>
                                            <td>MXNET-based deep neural networks up to 2 hidden layers, with multiple activation functions tested. Note that predefined grid search is used in hyperparameter optimization for this method (not random search). Details: <a href="https://github.com/kstawiski/OmicSelector/blob/master/R/OmicSelector_benchmark.R#L173" target="_blank">code</a>, <a href="https://mxnet.apache.org/" target="_blank">package</a>.
                                                <br>
                                                <i>It uses early stopping, but what maximum number of epochs should be used? <input class="form-control" id="search_iters_mxnet" name="search_iters_mxnet" type="number" min="1" max="500000" value="200" /></i>
                                            </td>
                                        </tr> -->
                                        <tr>
                                            <td><label class="switch"><input type="checkbox" name="method[]" value="xgbTree"><span class="slider round"></span></label></td>
                                            <td><code>xgbTree<code></td>
                                            <td>eXtreme gradient boosting. (note: this is a time-consuming method). Details: <a href="https://github.com/topepo/caret/blob/master/models/files/xgbTree.R" target="_blank">code</a>, <a href="https://www.rdocumentation.org/packages/xgboost" target="_blank">package</a>.</td>
                                        </tr>

                                        <tr>
                                            <td>Options:</td>
                                            <td colspan="2">
                                                <div class="form-group row">
                                                    <div class="col-sm-7">
                                                        <p><u>Number of random hyperparameter sets</u><br />
                                                            <font size="1"><i>(hyperparameter search is performed via <a href="https://topepo.github.io/caret/random-hyperparameter-search.html" target="_blank">random search</a>, how many should be checked?)</i></font>
                                                        </p>
                                                    </div>
                                                    <div class="col-sm-5">
                                                        <input class="form-control" id="search_iter" name="search_iter" type="number" min="1" max="500000" value="20" />
                                                    </div>
                                                </div>
                                                <div class="form-group row">
                                                    <div class="col-sm-7">
                                                        <p><u>Type of validation</u><br />
                                                            <font size="1"><i></i></font>
                                                        </p>
                                                    </div>
                                                    <div class="col-sm-5">
                                                        <select class="form-control" id="holdout" name="holdout">
                                                            <option value="TRUE" selected>Holdout-validation (on test set)</option>
                                                            <option value="FALSE">10-fold cross-validation repeated 5 times</option>
                                                        </select>
                                                    </div>
                                                </div>

                                            </td>
                                        </tr>
                                    </tbody>
                            </table>
                            <p>
                                <button type="submit" class="btn btn-success" value="Upload" name="submit" onclick="waitingDialog.show('Starting the analysis...');">
                                    <i class="fas fa-clipboard-check"></i>&emsp;Start benchmarking feature sets
                                </button>&emsp;</p>
                            </form>
                        </div>
                    </div>


                <?php
                    // Analiza gdy benchmark zrobiony.
                } else { ?>
                    <div class="panel panel-success">
                        <div class="panel-heading"><i class="fas fa-award"></i>&emsp;&emsp;Best signature</div>
                        <div class="panel-body">

                            <?php echo file_get_contents($target_dir . "best_signature.html"); ?>
                            <table class="table">
                                <tr>
                                    <td>Benchmark results:</td>
                                    <td>
                                        <p><a href="viewer.php?f=<?php echo $_GET['id']; ?>/benchmark.csv" class="btn btn-info" role="button" target="popup" onclick="window.open('viewer.php?f=<?php echo $_GET['id']; ?>/benchmark.csv','popup','width=1150,height=800'); return false;"><i class="fas fa-search-plus"></i> View</a>&emsp;<a href="/e/files/<?php echo $_GET['id']; ?>/benchmark.csv" class="btn btn-warning"><i class="fas fa-download"></i> Download</a>
                                            &emsp;
                                            <a href="viewer.php?f=<?php echo $_GET['id']; ?>/model_results.csv" class="btn btn-info" role="button" target="popup" onclick="window.open('viewer.php?f=<?php echo $_GET['id']; ?>/model_results.csv','popup','width=1150,height=800'); return false;"><i class="fas fa-search-plus"></i> View</a>&emsp;<a href="/e/files/<?php echo $_GET['id']; ?>/model_results.csv" class="btn btn-warning"><i class="fas fa-download"></i> Simplified</a>
                                            &emsp;<a href="/e/edit/<?php echo $_GET['id']; ?>/temp/benchmarkOmicSelector.log" class="btn btn-primary" role="button" target="popup" onclick="window.open('/e/edit/<?php echo $_GET['id']; ?>/temp/benchmarkOmicSelector.log','popup','width=1150,height=800'); return false;"><i class="fas fa-history"></i> View log</a>&emsp;<a href="/e/tree/<?php echo $_GET['id']; ?>/models" class="btn btn-primary" role="button" target="popup" onclick="window.open('/e/tree/<?php echo $_GET['id']; ?>/models','popup','width=1150,height=800'); return false;"><i class="fas fa-list"></i> Browse models</a></p>
                                    </td>
                                </tr>

                                <tr>
                                    <td>Customize this automatic analysis:</td>
                                    <td>
                                        <p><a href="/e/notebooks/<?php echo $_GET['id']; ?>/best_signature.Rmd" class="btn btn-danger" role="button" target="popup" onclick="window.open('/e/notebooks/<?php echo $_GET['id']; ?>/best_signature.Rmd','popup','width=1150,height=800'); return false;"><i class="fas fa-edit"></i> Edit & run</a>&emsp;<a href="process.php?type=best_signature_render&id=<?php echo $_GET['id']; ?>" class="btn btn-primary" role="button" onclick="waitingDialog.show('Rendering again.. this may take a while...');"><i class="fas fa-redo"></i> Refresh best signature analysis</a></p>
                                        <p><a href="process.php?type=delete_benchmark&analysisid=<?php echo $_GET['id']; ?>" class="btn btn-danger" role="button" onclick="return confirm('Are you sure? This will delete all the data regarding your benchmarking.')"><i class="fas fa-trash"></i> Re-run benchmarking (delete previous, model files will be intact)</a></p>
                                    </td>
                                </tr>
                            </table>

                        </div>
                    </div>
        </div>




    <?php } ?>
<?php } ?>

<div class="panel panel-warning">
    <div class="panel-heading"><i class="fas fa-puzzle-piece"></i>&emsp;&emsp;Post-analysis extensions</div>
    <div class="panel-body">
        <div class="panel panel-success">
            <div class="panel-heading clickable">
                <h3 class="panel-title" id="deep_learning">
                    <i class="fas fa-code-branch"></i>&emsp;<b>[DEEP LEARNING]</b> Training deep neural networks with grid search of hyperparameters.
                </h3>
            </div>
            <div class="panel-body">
                <p>Welcome to <b>deep learning</b> model training induction. This section will run deep learning extension and perform a grid search find the most optimal (deep) neural network for your datasets.</p>
                <form action="process.php?type=init_deeplearning" method="post" enctype="multipart/form-data"><input type="hidden" id="analysisid" name="analysisid" value="<?php echo $_GET['id']; ?>">
                    <div class="form-group">
                        <?php if (!file_exists($target_dir . "var_deeplearning_selected.txt")) { ?>
                            <p>Choose feature set the networks should be trained on <i>(note: only "all" features are avaiable if you do not run feature selection)</i>:
                                <select class="form-control" name="selected" id="selected">
                                    <option value="all">All features (without feature selection)</option>
                                    <?php
                                    if (file_exists($target_dir . "featureselection_formulas_final.csv")) {
                                        $types = array_map('str_getcsv', file($target_dir . "featureselection_formulas_final.csv"));
                                        $i = 1;
                                        foreach ($types as $row) {
                                            if ($i > 1) {
                                                echo  '<option value="' .  $row[0] . '">' .  $row[0] . ": " . $row[1] . '</option>';
                                            }
                                            $i = $i + 1;
                                        }
                                    } ?>
                                </select>
                            </p>
                            <p>Initial configuration of grid search <i>(note: autoencoders are expiermental)</i>:
                                <select class="form-control" name="autoencoders" id="autoencoders">
                                    <option value="0" <?php if (file('/PUBLIC', FILE_IGNORE_NEW_LINES)[0] == "1") { echo "disabled"; } ?>>Full scan: neural networks up to 3 hideen layers without autoencoders (97848 hyperparameter combinations)</option>
                                    <option value="1" <?php if (file('/PUBLIC', FILE_IGNORE_NEW_LINES)[0] == "1") { echo "disabled"; } ?>>Extendend scan: neural networks up to 3 hideen layers with and without autoencoders (293544 hyperparameter combinations)</option>
                                    <option value="2">Quick scan: neural networks with 1 hidden layers and without autoencoders (1994 hyperparameter combinations)</option>
                                </select>
                            </p>
                            <p>Use balanced training dataset <i>(if 'yes' the dataset balanced with SMOTE or ROSE will be used; the same balanced dataset as above is used)</i>:
                                <select class="form-control" name="balanced" id="balanced">
                                    <option value="FALSE">No. Use original dataset.</option>
                                    <option value="TRUE">Yes. Use balanced dataset.</option>
                                </select>
                            </p>
                        <?php } else { ?>
                            <p>You have already configured the deep learning training. If you wish to configure new initial settings click the button below. This is one-way step (i.e. you will not be able to resume training of new models in this set), but the results (models and settings) will be saved for the analysis.</p>
                            <p><a href="process.php?type=reconfigure_deep_learning&analysisid=<?php echo $_GET['id']; ?>" class="btn btn-primary" role="button"><i class="fas fa-people-arrows"></i> Save results and start new configuration</a></p>
                        <?php } 
                        if (file('/PUBLIC', FILE_IGNORE_NEW_LINES)[0] != "1") { ?>
                        <p>Number of parallel training processes <i>(this is highly dependent on your CPU/GPU and RAM, maximum should be picked based on trial and error, or just use e.g. 2 threads)</i>:
                            <input class="form-control" name="keras_threads" id="keras_threads" type="text" oninput="this.value=this.value.replace(/[^0-9]/g,'');" value="2" />
                        </p>
                        <?php } ?>

                        <button type="submit" class="btn btn-success" value="Upload" name="submit" onclick="waitingDialog.show('Setting up your enviorment...');">
                            <i class="fas fa-diagnoses"></i>&emsp;Start/resume deep learning
                        </button>
                        </div></form>

            </div>
        </div>


        <?php if (file_exists($target_dir . "deeplearning.csv")) { ?>
            <div class="panel panel-danger">
                <div class="panel-heading clickable">
                    <h3 class="panel-title" id="deep_learning_results">
                        <i class="fas fa-code-branch"></i>&emsp;<b>[DEEP LEARNING]</b> Deep neural networks - results.
                    </h3>
                </div>
                <div class="panel-body">
                    <table class="table">

                        <tbody>
                            <p>Files with results:</p>
                            <tr>
                                <td>Hyperparameters:</td>
                                <td><a href="viewer.php?f=<?php echo $_GET['id']; ?>/hyperparameters_deeplearning.csv" class="btn btn-info" role="button" target="popup" onclick="window.open('viewer.php?f=<?php echo $_GET['id']; ?>/hyperparameters_deeplearning.csv','popup','width=1150,height=800'); return false;"><i class="fas fa-search-plus"></i> View</a>&emsp;<a href="/e/files/<?php echo $_GET['id']; ?>/hyperparameters_deeplearning.csv" class="btn btn-warning"><i class="fas fa-download"></i> Download</a><br /><font size="1">Note: viewer may require significant amount of time to view it online. Downloading the file is recommended.</font></td>
                            </tr>
                            <tr>
                                <td>Performance of models (last configuration, main results):</td>
                                <td><a href="viewer.php?f=<?php echo $_GET['id']; ?>/deeplearning.csv" class="btn btn-info" role="button" target="popup" onclick="window.open('viewer.php?f=<?php echo $_GET['id']; ?>/deeplearning.csv','popup','width=1150,height=800'); return false;"><i class="fas fa-search-plus"></i> View</a>&emsp;<a href="/e/files/<?php echo $_GET['id']; ?>/deeplearning.csv" class="btn btn-warning"><i class="fas fa-download"></i> Download</a><br /><font size="1">Note: viewer may require significant amount of time to view it online. Downloading the file is recommended.</font></td>
                            </tr>
                            <tr>
                                <td>Merge:</td>
                                <td><a href="process.php?type=merge_deeplearning&analysisid=<?php echo $_GET['id']; ?>" class="btn btn-primary" role="button"><i class="fas fa-download"></i> Merge all deep learning runs and download results</a>
                                <br /><font size="1">Note: Merging is require for networks analysis and exporting.</font></td>
                            </tr>
                            <?php if (file_exists($target_dir . "merged_deeplearning.csv")) { ?>
                                <tr>
                                <td>Top 1000 networks (after merging):</td>
                                <td><a href="viewer.php?f=<?php echo $_GET['id']; ?>/merged_deeplearning_top.csv" class="btn btn-info" role="button" target="popup" onclick="window.open('viewer.php?f=<?php echo $_GET['id']; ?>/merged_deeplearning_top.csv','popup','width=1150,height=800'); return false;"><i class="fas fa-search-plus"></i> View</a>&emsp;<a href="/e/files/<?php echo $_GET['id']; ?>/merged_deeplearning_top.csv" class="btn btn-warning"><i class="fas fa-download"></i> Download</a>
                                <br /><font size="1">Note: Metaindex = mean of training, testing and validation accruacy. Metaindex2 = mean of testing and validaiton accuracy.</font></td>
                                </tr>
                                <tr>
                                <td>Choose best network:<br />
                                
                                </td>
                                <td>
                                    <p><font size=1">Choose network for further analysis:</font>
                                   
                                <div class="form-group">
                                <link href="https://cdn.jsdelivr.net/npm/select2@4.1.0-beta.1/dist/css/select2.min.css" rel="stylesheet" />
                                <script src="https://cdn.jsdelivr.net/npm/select2@4.1.0-beta.1/dist/js/select2.min.js"></script>
                                <script>$(document).ready(function () { $('#modelid').select2(); }); </script>
                                <form action="/deeplearning_model/" method="get" target="_blank">
                                <input type="hidden" id="analysisid" name="analysisid" value="<?php echo $_GET['id']; ?>">
                                <p><select class="form-control" name="modelid" id="modelid">
                                    <?php
                                    if (file_exists($target_dir . "merged_deeplearning_top.csv")) {
                                        $types = array_map('str_getcsv', file($target_dir . "merged_deeplearning_names.csv"));
                                        $i = 1;
                                        foreach ($types as $row) {
                                            if ($i > 1) {
                                                echo  '<option value="' .  $row[0] . '">' .  $row[0] . '</option>';
                                            }
                                            $i = $i + 1;
                                        }
                                    } ?>
                                </select></p>
                                <p><button type="submit" class="btn btn-success" value="submit" name="submit"><i class="fas fa-external-link-alt"></i>&emsp;Show neural network details and prediction tool</button></p>

                            </div></form></p>
                                </td>
                                </tr>


                            <?php } ?>

                        </tbody>
                    </table>

                </div>



            </div>
        <?php } ?>

    </div>
</div>

    <div class="panel panel-default">
        <div class="panel-heading"><i class="fas fa-bars"></i>&emsp;&emsp;Additional tools</div>
        <div class="panel-body"><a href="e/tree/<?php echo $_GET['id']; ?>" role="button" class="btn btn-primary" target="popup" onclick="window.open('/e/tree/<?php echo $_GET['id']; ?>','popup','width=1150,height=800'); return false;"><i class="fas fa-lock-open"></i>&emsp;Jupyter</a>&emsp;
            <a href="/process.php?type=rstudio&analysisid=<?php echo $_GET['id']; ?>" role="button" class="btn btn-primary" target="popup" onclick="window.open('/process.php?type=rstudio&analysisid=<?php echo $_GET['id']; ?>','popup','width=1150,height=800'); return false;"><i class="fas fa-lock-open"></i>&emsp;R Studio</a>&emsp;
            <a href="/process.php?type=radiant&analysisid=<?php echo $_GET['id']; ?>" role="button" class="btn btn-primary" target="popup" onclick="window.open('/process.php?type=radiant&analysisid=<?php echo $_GET['id']; ?>','popup3','width=1150,height=800'); return false;"><i class="fas fa-lock-open"></i>&emsp;Radiant</a>&emsp;
            <a href="/process.php?type=vscode&analysisid=<?php echo $_GET['id']; ?>" role="button" class="btn btn-primary" target="popup" onclick="window.open('/process.php?type=vscode&analysisid=<?php echo $_GET['id']; ?>','popup3','width=1150,height=800'); return false;"><i class="fas fa-lock-open"></i>&emsp;VS Code</a>&emsp;
            <a href="/" onclick="waitingDialog.show('Going back...');" class="btn btn-success"><i class="fas fa-sign-out-alt"></i>&emsp;Exit</a>
            <br><br><i>Login credentials to R Studio: username: <code><b><?php echo $_GET['id']; ?></b></code>, password: <code><b>OmicSelector</b></code></i>
        </div>
    </div>

    <!--Modal: Name-->
    <hr />
    <footer class="footer">
        <div class="container">
            <span class="text-muted">OmicSelector by Konrad Stawiski and Marcin Kaszkowiak&emsp;&emsp;&emsp;&emsp;<i class="fas fa-envelope"></i> konrad.stawiski@umed.lodz.pl&emsp;&emsp;&emsp;<i class="fas fa-globe-europe"></i>
                 <a href="https://biostat.umed.pl/OmicSelector/" target="_blank">https://biostat.umed.pl/OmicSelector/</a></span>
            <p>&emsp;</p>
        </div>
    </footer>
    <!-- /.container -->
    </div>
</body>

</html>