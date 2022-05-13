<?php
$fileName = "/OmicSelector/" . $_GET['f'];
if(!file_exists($fileName)) { die("File does not exists."); }


// Loading
if(!isset($_GET['loaded'])) {
  $actual_link = (isset($_SERVER['HTTPS']) && $_SERVER['HTTPS'] === 'on' ? "https" : "http") . "://$_SERVER[HTTP_HOST]$_SERVER[REQUEST_URI]";
  header('Location: /redirect.php?url='. urlencode($actual_link));
  die();
}

require 'vendor/autoload.php';

use PhpOffice\PhpSpreadsheet\Spreadsheet;
use PhpOffice\PhpSpreadsheet\Writer\Xlsx;

$reader = \PhpOffice\PhpSpreadsheet\IOFactory::createReaderForFile($fileName);
$reader->setReadDataOnly(true);
$spreadsheet = $reader->load($fileName);


//read excel data and store it into an array
$xls_data = $spreadsheet->getActiveSheet()->toArray(null, true, true, true);
/* $xls_data contains this array:
[1=>['A'=>'Domain', 'B'=>'Category', 'C'=>'Nr. Pages'], 2=>['A'=>'CoursesWeb.net', 'B'=>'Web Development', 'C'=>4000], 3=>['A'=>'MarPlo.net', 'B'=>'Courses & Games', 'C'=>15000]]
*/

//now it is created a html table with the excel file data
$html_tb ='<table class="table" id="table" name="table"><thead><tr><th>'. implode('</th><th>', $xls_data[1]) .'</th></tr></thead><tbody>';
$nr = count($xls_data); //number of rows
for($i=2; $i<=$nr; $i++){
  $html_tb .='<tr><td>'. implode('</td><td>', $xls_data[$i]) .'</td></tr>';
}
$html_tb .='</tbody></table>';
?>
<!DOCTYPE html>
    <html>
    <head>
    <title>OmicSelector - viewer</title>
    <link rel="stylesheet" href="https://maxcdn.bootstrapcdn.com/bootstrap/3.3.7/css/bootstrap.min.css"
      />
    <script src="https://code.jquery.com/jquery-3.4.1.min.js"  
    crossorigin="anonymous"></script>
    <link rel="stylesheet" href="https://maxcdn.bootstrapcdn.com/bootstrap/3.3.7/css/bootstrap-theme.min.css"
      crossorigin="anonymous" />
    <script src="https://maxcdn.bootstrapcdn.com/bootstrap/3.3.7/js/bootstrap.min.js"
     ></script>
    <meta charset="utf-8" />
    <meta http-equiv="X-UA-Compatible" content="IE=edge" />
    <meta name="viewport" content="width=device-width, initial-scale=1" />
    <meta name="description" content="" />
    <meta name="author" content="" />
    <script src="https://cdnjs.cloudflare.com/ajax/libs/font-awesome/5.12.1/js/all.min.js"  ></script>
    <script src="https://code.jquery.com/jquery-1.8.3.min.js" type="text/javascript"></script>
    <script src="https://code.jquery.com/ui/1.9.2/jquery-ui.js" type="text/javascript"></script>
        
    <script src="https://cdn.datatables.net/1.10.20/js/jquery.dataTables.min.js" type="text/javascript"></script>
    <link rel="stylesheet" href="https://cdn.datatables.net/1.10.20/css/jquery.dataTables.min.css" crossorigin="anonymous" /> 
        <script src="https://code.jquery.com/jquery-3.3.1.js
<script src="https://cdn.datatables.net/1.10.20/js/jquery.dataTables.min.js" type="text/javascript"></script>
<script src="https://cdn.datatables.net/buttons/1.6.1/js/dataTables.buttons.min.js" type="text/javascript"></script>
<script src="https://cdn.datatables.net/buttons/1.6.1/js/buttons.flash.min.js" type="text/javascript"></script>
<script src="https://cdnjs.cloudflare.com/ajax/libs/jszip/3.1.3/jszip.min.js" type="text/javascript"></script>
<script src="https://cdnjs.cloudflare.com/ajax/libs/pdfmake/0.1.53/pdfmake.min.js" type="text/javascript"></script>
<script src="https://cdnjs.cloudflare.com/ajax/libs/pdfmake/0.1.53/vfs_fonts.js" type="text/javascript"></script>
<script src="https://cdn.datatables.net/buttons/1.6.1/js/buttons.html5.min.js" type="text/javascript"></script>
<script src="https://cdn.datatables.net/buttons/1.6.1/js/buttons.print.min.js" type="text/javascript"></script>
        <script>
        $(document).ready( function () {
    $('#table').DataTable({
    dom: 'Bfrtip',
        buttons: [
            'copy', 'csv', 'excel', 'pdf', 'print'
        ],
        "pageLength": 50
    });
} );
        </script>
    </head>
    <body>
<?php echo $html_tb; ?>
    </body>
    </html>