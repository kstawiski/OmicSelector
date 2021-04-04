<html>
  <head>
    <title>OmicSelector - system monitor</title>
    <link rel="stylesheet" href="https://maxcdn.bootstrapcdn.com/bootstrap/3.3.7/css/bootstrap.min.css"
      />
    <meta http-equiv="refresh" content="5" />
  </head>
  <body>
		<pre><?php passthru('/usr/bin/nvidia-smi'); ?></pre>
  </body>
</html>