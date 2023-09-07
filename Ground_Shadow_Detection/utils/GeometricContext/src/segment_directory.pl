#!/usr/bin/perl -w



####################

$input_dir = 'C:/Users/hasanbalci/Documents/MATLAB/Ground_Shadow_Detection/GeometricContext/src//';
$output_dir = 'C:/Users/hasanbalci/Documents/MATLAB/Ground_Shadow_Detection/GeometricContext/src/';

$convertfrom = '.ppm';

@imagelist = glob $input_dir.'*'.$convertfrom;
print $input_dir.'*'.$convertfrom . "\n";

$count = 0;
foreach $image (@imagelist) {    
    $count = $count + 1;
    if ($count > 0) {
	$where = index($image, $convertfrom); 
	$where2 = index($image, '//');    
	$basename = substr($image, $where2+2, $where-$where2-2);
   
	print $basename .  "\n";
     	system('segment 0.8 100 100 ' . $image . ' ' . $output_dir . $basename . '.ppm');
    }

}
