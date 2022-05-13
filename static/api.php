<?php
if(!isset($_GET['typ'])) { die("Nie podano typu"); }

function tailCustom($filepath, $lines = 1, $adaptive = true) {

    // Open file
    $f = @fopen($filepath, "rb");
    if ($f === false) return false;

    // Sets buffer size, according to the number of lines to retrieve.
    // This gives a performance boost when reading a few lines from the file.
    if (!$adaptive) $buffer = 4096;
    else $buffer = ($lines < 2 ? 64 : ($lines < 10 ? 512 : 4096));

    // Jump to last character
    fseek($f, -1, SEEK_END);

    // Read it and adjust line number if necessary
    // (Otherwise the result would be wrong if file doesn't end with a blank line)
    if (fread($f, 1) != "\n") $lines -= 1;
    
    // Start reading
    $output = '';
    $chunk = '';

    // While we would like more
    while (ftell($f) > 0 && $lines >= 0) {

        // Figure out how far back we should jump
        $seek = min(ftell($f), $buffer);

        // Do the jump (backwards, relative to where we are)
        fseek($f, -$seek, SEEK_CUR);

        // Read a chunk and prepend it to our output
        $output = ($chunk = fread($f, $seek)) . $output;

        // Jump back to where we started reading
        fseek($f, -mb_strlen($chunk, '8bit'), SEEK_CUR);

        // Decrease our line counter
        $lines -= substr_count($chunk, "\n");

    }

    // While we have too many lines
    // (Because of buffer size we might have read too many)
    while ($lines++ < 0) {

        // Find first newline and remove all text before that
        $output = substr($output, strpos($output, "\n") + 1);

    }

    // Close file and return
    fclose($f);
    return trim($output);

}

switch ($_GET['typ']) {
    case "encrypt":
        $ciphering = "BF-CBC"; 
        $iv_length = openssl_cipher_iv_length($ciphering); 
        $options = 0; 
        $encryption_iv = random_bytes($iv_length); 
        $encryption_key = openssl_digest(php_uname(), 'MD5', TRUE); 
        $decryption_iv = random_bytes($iv_length); 
        $decryption_key = openssl_digest(php_uname(), 'MD5', TRUE); 
        echo openssl_encrypt($_GET['txt'], $ciphering, $encryption_key, $options, $encryption_iv); 
        break;
    case "decrypt":
        $ciphering = "BF-CBC"; 
        $iv_length = openssl_cipher_iv_length($ciphering); 
        $options = 0; 
        $encryption_iv = random_bytes($iv_length); 
        $encryption_key = openssl_digest(php_uname(), 'MD5', TRUE); 
        $decryption_iv = random_bytes($iv_length); 
        $decryption_key = openssl_digest(php_uname(), 'MD5', TRUE); 
        echo openssl_decrypt($_GET['txt'], $ciphering, $decryption_key, $options, $encryption_iv); 
        break;
    case "odczyt":
        if(file_exists($_GET['plik'])) {
        //echo file_get_contents($_GET['plik']);
        echo tailCustom($_GET['plik'], 1000);
        } else { echo "Incorrect request"; }
        break;
    case "cpu":
        passthru('/usr/bin/top -b -n 1');
        break;
    case "gpu":
        // passthru('/usr/bin/nvidia-smi');
        passthru('/bin/bash -c \'source ~/.bashrc && nvidia-smi\'');
        break;

}



die();
?>