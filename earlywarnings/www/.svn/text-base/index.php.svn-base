
<!-- This is the project specific website template -->
<!-- It can be changed as liked or replaced by other content -->

<?php

$domain=ereg_replace('[^\.]*\.(.*)$','\1',$_SERVER['HTTP_HOST']);
$group_name=ereg_replace('([^\.]*)\..*$','\1',$_SERVER['HTTP_HOST']);
$themeroot='r-forge.r-project.org/themes/rforge/';

echo '<?xml version="1.0" encoding="UTF-8"?>';
?>
<!DOCTYPE html
	PUBLIC "-//W3C//DTD XHTML 1.0 Transitional//EN"
	"http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd">
<html xmlns="http://www.w3.org/1999/xhtml" xml:lang="en" lang="en   ">

  <head>
	<meta http-equiv="Content-Type" content="text/html; charset=UTF-8" />
	<title><?php echo $group_name; ?></title>
	<link href="http://<?php echo $themeroot; ?>styles/estilo1.css" rel="stylesheet" type="text/css" />
  </head>

<body>

<!-- R-Forge Logo -->
<table border="0" width="100%" cellspacing="0" cellpadding="0">
<tr><td>
<a href="http://r-forge.r-project.org/"><img src="http://<?php echo $themeroot; ?>/imagesrf/logo.png" border="0" alt="R-Forge Logo" /> </a> </td> </tr>
</table>


<!-- get project title  -->
<!-- own website starts here, the following may be changed as you like -->

<?php if ($handle=fopen('http://'.$domain.'/export/projtitl.php?group_name='.$group_name,'r')){
$contents = '';
while (!feof($handle)) {
	$contents .= fread($handle, 8192);
}
fclose($handle);
echo $contents; } ?>

<!-- end of project description -->

<p>Critical transitions are abrupt sudden changes in the dynamics of systems. Recent theoretical insights suggest that specific statistical properties may change in predictable manners prior to these critical transitions. The Early Warning Signals Toobox summarizes available methods and approaches for identifying critical transitions in timeseries data.</p>
<p>&nbsp;</p>
<blockquote>
  <blockquote>
    <blockquote>
      <blockquote>
        <blockquote>
          <p><img src="images/fig_website_R.jpg" alt="" name="figure" width="572" height="501" align="top" id="figure" /></p>
        </blockquote>
      </blockquote>
    </blockquote>
  </blockquote>
</blockquote>
<p>&nbsp;</p>
<p>The toolbox is currently under constant development. We aim in including all available methods and indicators that have been proposed as potential tools for signaling upcoming critical transitions.</p>
<p><strong>References</strong></p>
<p><a name="_ENREF_26" id="_ENREF_26">Scheffer, M., J. Bascompte, W. A. Brock, V. Brovkin, S. R. Carpenter, V. Dakos, H. Held, E. H. van Nes, M. Rietkerk, and G. Sugihara. 2009. Early-warning signals for critical transitions. Nature 461:53-59.</a></p>
<p><a name="_ENREF_12" id="_ENREF_12">Dakos, V., M. Scheffer, E. H. van Nes, V. Brovkin, V. Petoukhov, and H. Held. 2008. Slowing down as an early warning signal for abrupt climate change. Proceedings of the National Academy of Sciences 105:14308-14312.</a></p>
<p> The <strong>project summary page</strong> you can find <a href="http://<?php echo $domain; ?>/projects/<?php echo $group_name; ?>/"><strong>here</strong></a>. </p>

</body>
</html>
