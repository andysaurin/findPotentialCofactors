#!/usr/bin/php -n
<?php
/**********
This script uses MOODS to find PWM occurrences:
http://www.cs.helsinki.fi/group/pssmfind/

J. Korhonen, P. Martinmäki, C. Pizzi, P. Rastas and E. Ukkonen. MOODS: fast search for position weight matrix matches in DNA sequences. Bioinformatics 25(23), pages 3181-3182. (2009)

forked at https://github.com/andysaurin/MOODS
**********/

require_once( dirname(__FILE__) . "/args.php" );
require_once( dirname(__FILE__) . "/bioPHP/reader_gff_fasta.php" );
require_once( dirname(__FILE__) . "/bioPHP/seq_alignment.php" );

// path to find_pssm_dna
// https://github.com/andysaurin/MOODS
$find_pssm_dna = dirname(__FILE__).'/find_pssm_dna';

$uniqid = uniqid();
$defaultP = '0.001';
$default_distance = null;
$output_file = 'output.txt';
$align_sequences = false;
$plot = true;

$help_msg = "\nUsage:\n{$args[exec]} --seq=[fasta_file] --m1==[matrix1.pwm] --p1=[p-value] --m2==[matrix2.pwm] --p2=[p-value] --d=100

	Options:
	--seq		Required : The input FASTA sequence file to compare the motifs against
	--m1		Required : First PWM matrix to test with in TAB format
	--p1		Optional : The minimum p-value for scoring the First PWM against the sequence for it to be considered a hit location (Default: {$defaultP})
	--m2		Required : Second PWM matrix to test with in TAB format
	--p2		Optional : The minimum p-value for scoring the Second PWM against the sequence for it to be considered a hit location (Default: {$defaultP})
	--d		Optional : The maximum distance (in base pairs) between PWM hit locations to be considered a PWM-pair match (Default: No minimum)
	--o		Optional : The output filename (Default: {$output_file})
	--aln	Optional : Perform sequence alignment on matches (Default: false)
	--plot	Optional : Report hits as distance hits to use in scatter plot (Default : false)

	Note:
		1) Note the PWMs must be in tab-format. Use RSAT to convert

";


if ( !is_executable($find_pssm_dna) )
	die("{$find_pssm_dna} is not executable.\n");

if ( !$argv || count($argv) < 2 )
	die($help_msg);




/*****
Options test
*****/
//sequence file
if ( isset($args['options']['seq']) ) {
	$fasta = $args['options']['seq'];

	if ( !is_file($fasta) || !is_readable($fasta) )
		$errors[] = "\tThe FASTA sequence file {$fasta} down not exist or cannot be read.\n";

} else {
	$errors[] = "\tThe FASTA sequence file is not set: --seq=<fasta_file>\n";
}

//p-values
if ( isset($args['options']['p1']) ) {
	$pValues['m1'] = $args['options']['p1'];
} else {
	$pValues['m1'] = $defaultP;
}
if ( isset($args['options']['p2']) ) {
	$pValues['m2'] = $args['options']['p2'];
} else {
	$pValues['m2'] = $defaultP;
}

if ( isset($args['options']['d']) ) {
	$min_dist = $args['options']['d'];
} else {
	$min_dist = $default_distance;
}

//output file
if ( isset($args['options']['o']) ) {
	$output_file = $args['options']['o'];
}


if ( $pValues['m1'] >=1 || $pValues['m1'] <= 0 || $pValues['m2'] >=1 || $pValues['m2'] <= 0)
	$errors[] = "\tThe p-values must be greater than 0 and less than 1.\n";


if ( !isset($args['options']['m1']) || !isset($args['options']['m2']) )
	$errors[] = "\tYou must supply at least two PWM motifs files\n";
else {
	$pwms['m1'] = $args['options']['m1'];
	$pwms['m2'] = $args['options']['m2'];


	foreach ($pwms as $k=> $pwm) {
		if ( !is_file($pwm) || !is_readable($pwm) ) {
			$errors[] = "\tThe PWM matrix file {$pwm} down not exist or cannot be read.\n";
		} else {
			//check format
			$tmp_pwm = file_get_contents($pwm);
			$tmp_pwm = str_replace("[", "", $tmp_pwm);
			$tmp_pwm = str_replace("]", "", $tmp_pwm);
			$tmp_pwm = preg_replace("/[ ]/", "\t", $tmp_pwm); //spaces to tabs
/*
			$tmp_pwm = preg_replace("/\n[\t]+/", "\n", $tmp_pwm); //remove starting tabs
			$tmp_pwm = preg_replace("/^[\t]+/", "", $tmp_pwm); //remove starting tabs
*/
			$tmp_pwm = preg_replace("/[^0-9\t\n\.]/i", "", $tmp_pwm);

			$matrix = explode("\n", $tmp_pwm);
			if ( is_array($matrix) && count($matrix) ) { //trim the matrix file of whitespace
				foreach($matrix as $k2=>$v) {
					$matrix[$k2] = trim($v);
				}
				$tmp_pwm = implode("\n", $matrix);

			} else {
				$matrix = array();
			}
			preg_match_all("|\n|", $tmp_pwm, $matches);

//print_r($matches);
			if ( !preg_match("/^[0-9\t\n\.]+$/", $tmp_pwm) || count($matches[0]) != 4 ) {

				$errors[] = "\tThe PWM matrix file {$pwm} is in the wrong format. PWMs must me in TAB format eg:
0       0       0       173     0       0       173     57
173     0       0       0       0       173     0       32
0       0       0       0       0       0       0       62
0       173     173     0       173     0       0       22

Given PWM is:

{$tmp_pwm}
\n";
print_r($matches);
			} else {

				$tmp = sys_get_temp_dir() . '/pwm-' . mt_rand() . '-' .md5($pwm);
				$fh = fopen($tmp, "w");

				fwrite($fh, $tmp_pwm);
				$tmp_pwms[$k]['tmp_loc'] = $tmp;
				$tmp_pwms[$k]['pValue'] = $pValues[$k];
				$tmp_pwms[$k]['name'] = pathinfo($pwm, PATHINFO_FILENAME);
				//get the length
				$pwm_array = explode("\n", $tmp_pwm);

				preg_match_all("|\t+|", $pwm_array[0], $matches);

				$tmp_pwms[$k]['length'] = count($matches[0]) + 1;

				//reverse the pwm
				$reverse_pwm = reverse_complement_pwm($tmp_pwm);

				$fh = fopen($tmp."-rc", "w");
				fwrite($fh, $reverse_pwm);

//$tmp_pwms[$k]['length'] = 0;
			}
		}
	}
}

//print_r($tmp_pwms);

if ( isset($errors) && is_array($errors) && count($errors) ) {

	fwrite(STDERR, $help_msg."The following errors were encountered:\n");
	foreach ($errors as $error) {
		fwrite(STDERR, $error);
	}
	fwrite(STDERR, "\n");
	exit(1);
}


$fa = new Reader();
$fa->INPUT_FILE = $fasta;
$fa->TYPE_FILE = 'fasta';
$fa->read();

/**
 when finding the PWMs all the sequences are nitted together to give one long sequence

 Thus, we need to set up an array of the start positions of each sequence
**/

$total_bases = 0;
//$seq_ids = array();

$next_key = 1;
/*
print_r($fa->OBJECT[0]);
exit;
*/

$longest_sequence = 0;
$cat_sequences = ''; //this will be all the sequences catenated together


//print_r($fa->UID);

foreach ( $fa->OBJECT as $k=> $fa_object ) {

//echo strlen($fa_object->getSequence());
//print_r($fa_object);
	$seq = 	$fa_object->getSequence();
	$len = strlen($seq) - 1; //the seq is always less 1 (new line char)

	if ($len > $longest_sequence)
		$longest_sequence = $len; //this is used when scanning distance between motifs when --d is not set

	$total_bases = $total_bases + $len;

	$seq_ids[$next_key] = $fa_object->getId();

	$next_key = ($next_key + $len); //set the next key

	$cat_sequences .= trim($seq);
}

/*
if ($min_dist === null)
	$min_dist = $longest_sequence;
*/
/*
print_r($seq_ids);
gives:
Array
(
    [1] => enriched_region_93
    [2283] => enriched_region_200
    [2905] => enriched_region_201
    [6939] => enriched_region_246
    [8447] => enriched_region_284
...
)

*/

//free up memory
//$fa = '';

//now do the PWM search on the FASTA file

// find_pssm_dna 0.001 expressed.fa srp.pwm

fwrite(STDERR, "+++++++++++++++++\nLoaded ".count($seq_ids)." sequences totalling {$total_bases} bases.\n+++++++++++++++++\n");

foreach ($tmp_pwms as $k=> $pwm ) {

	$pValue = $pwm['pValue'];

	$command = "{$find_pssm_dna} {$pValue} {$fasta} {$pwm[tmp_loc]} {$pwm[tmp_loc]}-rc";
	$matches = `{$command} 2>/dev/null`;
//echo $command."\n";

	$matches = trim($matches);

	$hits = explode("\n", $matches);

	fwrite(STDERR, "{$pwm[name]} yielded ".count($hits)." hits having a p-value of " . $pwm['pValue'] . " or less.\n");
	@unlink($pwm['tmp_loc']);
	@unlink($pwm['tmp_loc']."-rc");

	if ( count($hits) > 0 ) {

		foreach ($hits as $hk =>$hit) {
			$array = explode("\t", $hit); // $hit is <hitPosition>	<score>

			$position = $array[0]; //the base position in the flattened sequence file
			$tmp_pwms[$k]['hits'][$position] = $array[1]; //the score
/*
			if ($array[1] < 8) {
				unset($hits[$hk]);
			} else {

				$position = $array[0]; //the base position in the flattened sequence file
				$tmp_pwms[$k]['hits'][$position] = $array[1]; //the score
			}
*/

		}

	}else
		$no_hits = 1;


}

if ( $no_hits == 1 ) {
	fwrite(STDERR, "At least one of the PWMs yielded no hits. Terminating.\n");
	exit();
}
//now we have the PWM hits (and their scores), lets see if any hits are close by

if ($min_dist !== null )
	fwrite(STDERR, "\nScanning for sequences containing both PWMs closer than {$min_dist}bp apart\n");
else
	fwrite(STDERR, "\nScanning for sequences containing both PWMs\n");

foreach ($tmp_pwms as $k=> $pwm ) {

	$pwm_name = $pwm['name'];
	if ( !is_array($pwm['hits']) )
		break; //PWM yielded no hits

	//group the hits into the DNA regions
	foreach($pwm['hits'] as $pos=>$score) {

		//find the associated region_id by the starting position of the hit and counting backwards till you hit the region starting base
		for($i=$pos; $i>0; $i--) {
			if ( isset($seq_ids[$i]) ) { //yup, we hit it
				$region_name = $seq_ids[$i];
				$regions[$region_name][$pwm_name][$pos] = $score;
				break;
			}
		}

	}


}

if ( !is_array($regions) )
	die("No PWM matches were found in the sequence\n"); //should never equate, but safe than sorry


foreach ($regions as $region_name=> $array) {


/*
print_r($array);
Array
(
    [abda_dmmpmm] => Array
        (
            [1105] => 4.40672
            [1113] => 4.40672
        )

    [srp] => Array
        (
            [2046] => 3.48614
        )

)
*/

	foreach($tmp_pwms as $key=>$pwm){
		$pwm_name = $pwm['name'];

		if ( !is_array($array[$pwm_name]) ) { //only one PWM in this region - so it's no use
			unset($regions[$region_name]);

 			break;
		}

	}

}


//this is the list of regions with at least 1 of each PWM in it
foreach ($regions as $region_name=> $array) {

	$i = 0;

	$matches = null;
	$tmp_array = array();
	$hits = null;

	foreach($tmp_pwms as $key=>$pwm){
		$pwm_name = $pwm['name'];
		$tmp_array[] = $array[$pwm_name]; //this screws up the array naming, but the problem was naming the key as the PWM name
	}

	$matches = find_closest($tmp_array);

	if (is_array($matches) ) {


		foreach ($matches as $pos1=> $tmp_array) {

			$hits[$region_name][$i]['distance'] = $tmp_array['closest_distance'];
			$hits[$region_name][$i]['score'] = $tmp_array['score'];

			$pos2 = $tmp_array['closest_position'];

			//get the PWM name of the 1st position
			foreach($tmp_pwms as $key=>$pwm){

				$pwm_name = $pwm['name'];

				if ( isset($array[$pwm_name][$pos1]) ) {
					$hits[$region_name][$i][$pwm_name] = $pos1;
				}
				if ( isset($array[$pwm_name][$pos2]) ) {
					$hits[$region_name][$i][$pwm_name] = $pos2;
				}
			}

			$i++;
		}
		//$hits[$region_name][$i] =
//		echo "yeah!\n";



//		print_r($matches);
//		print_r($array);

		$hits[$region_name] = unique_pwm_pairs($hits[$region_name]);

		$total_hits = $total_hits + count($hits[$region_name]);
		$total_hits_array[] = $hits;

//		exit();
	}

}

//print_r($total_hits_array);

if ( $total_hits < 1 )
	die("No DNA sequences were found.\n");

if ($min_dist != null && $min_dist != $longest_sequence )
	fwrite(STDERR, "\n{$total_hits} DNA regions out of ".count($seq_ids)." contain PWMs within {$min_dist}bp of each other.\n");
else
	fwrite(STDERR, "\n{$total_hits} DNA regions out of ".count($seq_ids)." contain at least one of each PWMs.\n");

if ( $align_sequences === true )
	fwrite(STDERR, "\nPerforming sequence alignments (this may take some time)...\n");

foreach ($total_hits_array as $k=>$region_array) {

	foreach ($region_array as $region_name=>$array1) {
		foreach($array1 as $k1=>$array2) {

			if ( $align_sequences === true ) {
				$total_hits_array[$k][$region_name][$k1] = get_sequence($array2, $region_name);
			}

			$distances[$k] = $array2['distance'];
			$scores[$k] = log(( (1/$array2['distance']) * $array2['score'] ), 2);
		}
	}

}
//$total_hits_array = get_sequence($total_hits_array);

$fh = fopen($output_file, "w+");

if ($plot === true) {

	$max_dist = max($distances);
	$min_dist = min($distances);

	$frequencies = array_count_values($distances);

/*
	$output = "Distance	Count\n";
	for ($x=$min_dist; $x<=$max_dist; $x++) {
		$output .= "{$x}\t".$distances[$x]."\n";
	}
*/
/*
	sort($distances);
	sort($scores);
*/

/*
	$output = implode(",", $distances);
*/
/*
	$output = implode(",", $scores);
	$output = trim($output);
*/

	foreach ($distances as $k=>$distance) {
		$output .="{$distance}\t{$scores[$k]}\n";
	}
} else {

	$output = var_export($total_hits_array, true);

}
fwrite($fh, $output);

fwrite(STDERR, "\nAnalysis Complete.\n\tOutput has been written to {$output_file}\n\n");
exit();

function get_sequence($array, $region_name=null){

	if ( !is_array($array) )
		return null;

	if ( !isset($array['distance']) )
		return get_sequence($array);

	global $cat_sequences;
	global $tmp_pwms;
	global $fa;
	global $align_sequences;

	//find which is the start sequence
	$end_pos = 0;
	foreach ($tmp_pwms as $pwm) {
		$pwm_name = $pwm['name'];


		if ( !isset($start_pos) || $array[$pwm_name] < $start_pos) {
			$start_pos = $array[$pwm_name];
		}

		if (  ($array[$pwm_name] + $pwm['length']) > $end_pos ) {
			$end_pos = $array[$pwm_name] + $pwm['length'];
		}

	}

	$start_pos = $start_pos;
	$offset = $end_pos - $start_pos;
	//$end_pos = $start_pos + $array['distance'];

	$array['sequence'] = substr($cat_sequences, $start_pos, $offset );

	if ( $align_sequences === false )
		return $array;

	if ( !$region_name) {
		$array['alignment'] = null;
		return $array;
	}
	$uid = array_search($region_name, $fa->UID);
	if (!is_numeric($uid)) {
		$array['alignment'] = null;
		return $array;
	}

	//perform sequence alignment
	$region_sequnce = $fa->OBJECT[$uid]->getSequence();

	$alignment=align_DNA($array['sequence'],$region_sequnce);

	$align_seqa=$alignment["seqa"];
	$align_seqb=$alignment["seqb"];

	// COMPARE ALIGNMENTS
	$compare=compare_alignment($align_seqa,$align_seqb);

	$aln = "\n";
	$i=0;
	while($i<strlen($align_seqa)){
	        $ii=$i+100;
	        if ($ii>strlen($align_seqa)){$ii=strlen($align_seqa);}
	        $aln .= substr($align_seqa,$i,100)."  $ii\n";
	        $aln .= substr($compare,$i,100)."\n";
	        $aln .= substr($align_seqb,$i,100)."  $ii\n\n";
	        $i+=100;
	}


	$array['alignment'] = $aln;

	return $array;

}


function unique_pwm_pairs($hits) {

	if ( count($hits) == 1 )
		return $hits; //nothing to do here as there's only one PWM pair


	global $tmp_pwms;

/**

We need to remove duplicate entries - ie below [abda_dmmpmm] => 453182 has been matched twice, so keep only the shortest distance

[enriched_region_9517] => Array
    (
        [0] => Array
            (
                [distance] => 157
                [abda_dmmpmm] => 453182
                [srp] => 453339
            )

        [1] => Array
            (
                [distance] => 186
                [abda_dmmpmm] => 453182
                [srp] => 453368
            )

    )

**/
//print_r($hits);
	foreach ($hits as $hk=>$hv) {

		foreach($tmp_pwms as $key=>$pwm){

			$pwm_name = $pwm['name'];
			$pwm_pos = $hv[$pwm_name];
			$pwm_distance = $hv['distance'];

			foreach ($hits as $hk2=>$hv2) {
				if ( $hk2 != $hk ) { //don't check ourself from the outer foreach

					if ( $hv2[$pwm_name] == $pwm_pos ) {

						//are we a longer distance than the the outer foreach pwm?
						if ( $hv2['distance'] > $pwm_distance ) {

							unset($hits[$hk2]);

							if (count($hits) == 1)
								return $hits;
							else
								continue;

						} else {

							unset($hits[$hk]);
							break(2);
/*
print_r($hits);
print_r($hv2);
							die('here2');
*/
						}

					}

				}

			}


		}
	}

	return $hits;

}


function find_closest($array1, $array2=null) {

	global $min_dist;
	global $longest_sequence;

	if ( !$array2 ) {
//die('here');
		return find_closest($array1[0], $array1[1]);

	}
//die( print_r($array2) );

	if ( count($array1) > count($array2) ) {
		//wrong way around - we want to check the smallest array against the largest
		return find_closest($array2, $array1);
	}

	//$i = 0;
	foreach ($array1 as $pos=>$score) {

		$closest_distance = null;
		$closest_position = null;

		$score_pwm1 = $score; //the score for this matrix hit

		//find nearest hit in $array2
		if ($min_dist !== null ) {

			$start_pos = $pos - $min_dist;
			$end_pos = $pos + $min_dist;

		} else { //as we have no idea of the DNA size, we'll use the longest piece of DNA to itereate through

			$start_pos = $pos - $longest_sequence;
			$end_pos = $pos + $longest_sequence;

		}

		for ($i=$start_pos; $i<=$end_pos; $i++) {

			if ( isset($array2[$i]) ) {
				$hit_pos = $i;

				$score_pwm2 = $array2[$i]; //the score for the partner matrix hit

				if ($hit_pos > $pos)
					$hit_distance = $hit_pos - $pos;
				else
					$hit_distance = $pos - $hit_pos;

				if ( !$closest_distance ) {
					$closest_distance = $hit_distance;
					$closest_position = $hit_pos;
				} else {
					if ( $hit_distance < $closest_distance ) {
						$closest_distance = $hit_distance;
						$closest_position = $hit_pos;
					}

				}
			}
		}

		if ( $min_dist == null || ($closest_distance && $closest_distance <= $min_dist) ) { //yeah, we have 2 PWMs closer than the max distance spacing

			$hits[$pos]['closest_position'] =  $closest_position;
			$hits[$pos]['closest_distance'] =  $closest_distance;
			$hits[$pos]['score'] = $score_pwm1 + $score_pwm2; //the combined score

		}

	}
	if ( isset($hits) && is_array($hits) ) {
		return $hits;
	}else
		return false;

}

//print_r($regions);
//print_r($total_hits_array);

function reverse_complement_pwm($pwm) {

	$array = explode("\n", $pwm);

	$array = array_reverse($array);

	foreach ($array as $k=>$v) {
		$pieces = explode("\t", $v);
		$pieces = array_reverse($pieces);
		$array[$k] = implode("\t", $pieces);
	}

	$pwm_rc = implode("\n", $array);

	return $pwm_rc;

}
?>