<html>
  <head>
    <title>OmicSelector - system monitor</title>
    <link rel="stylesheet" href="https://maxcdn.bootstrapcdn.com/bootstrap/3.3.7/css/bootstrap.min.css"
    integrity="sha384-BVYiiSIFeK1dGmJRAkycuHAHRg32OmUcww7on3RYdg4Va+PmSTsz/K68vbdEjh4u" crossorigin="anonymous" />
    <meta http-equiv="refresh" content="2" />
  </head>
  <body>
		<pre><?php passthru('/usr/bin/top -b -n 1'); ?></pre>
  </body>
</html>