<html>
  <head>
    <title>OmicSelector - CPU monitor</title>
    <link rel="stylesheet" href="https://maxcdn.bootstrapcdn.com/bootstrap/3.3.7/css/bootstrap.min.css"
      />
    <script
  src="https://code.jquery.com/jquery-3.6.0.js"
  integrity="sha256-H+K7U5CnXl1h5ywQfKtSj8PCmoN9aaq30gDh27Xc0jk="
  crossorigin="anonymous"></script>
  </head>
  <body>
  <script>
    $( document ).ready(function() {

    function get_log(){
        var feedback = $.ajax({
            type: "GET",
            url: "api.php",
            async: false,
            data: {
                "typ" : "cpu"
            },
            success: function() {
                //setTimeout(function(){ get_log();}, 1000);
            }
        }).responseText;
        
        $('#log').html(feedback);
        $('#log').animate({scrollTop:document.getElementById("log").scrollHeight}, 'slow');
    }
    
    get_log();
    setInterval(get_log, 2000);
    
    });
</script>
		<pre id="log"></pre>
  </body>
</html>