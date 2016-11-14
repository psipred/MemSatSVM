#!/usr/bin/perl
## MEMSAT-SVM script
##
## *********************************************************
## *    MEMSAT-SVM - MEMbrane protein Structure            *
## *    And Topology using Support Vector Machines         *
## * Integral Membrane Protein Topology Prediction Program *
## *    Copyright (C) 2008 Tim Nugent                      *
## *********************************************************
##
## This program is copyright and may not be distributed without
## permission of the author unless specifically permitted under
## the terms of the license agreement.
##
## THIS SOFTWARE MAY ONLY BE USED FOR NON-COMMERCIAL PURPOSES. PLEASE CONTACT
## THE AUTHOR IF YOU REQUIRE A LICENSE FOR COMMERCIAL USE.

use strict;
use warnings;
use Cwd;
use Getopt::Long;

## IMPORTANT : this variable must be set here to the absolute path to the
## library called "lib" and found in the same directory as this script.
use FindBin;
use lib "$FindBin::Bin/lib";
use Math::Round qw(:all);

## IMPORTANT : this variable must be set to the directory containing this
## script either here or using flag '-w' at runtime.
# my $mem_dir = '/cs/research/bioinf/home1/green/dbuchan/Code/memsat-svm';
my $mem_dir = '';

## IMPORTANT : these paths (folder with the NCBI executables and path to the
## database) must be set either here or using the appropriate flags at runtime.
# my $ncbidir = '/scratch0/NOT_BACKED_UP/dbuchan/Applications/blast-2.2.26/bin/';
# my $dbname = '/scratch0/NOT_BACKED_UP/dbuchan/uniref/test_db.fasta';
my $ncbidir = '';
my $dbname = '';

## Executable and input/output paths are initialised in sub get_arguments.
my ($input_path, $output_path, $model_path, $datadir, $svm_classify,
    $memsat_svm_bin_path, $globmem_bin_path, $mem_pred_svm_bin_path, $nnsat_bin_path);
my $mtx = 0;
my $remove_files = 0;
my $globmem = 0;
my $graphics = 1;
my $format = 1;
my $runmem3 = 0;
my $cores = 1;
my $globmem_score_threshold = -0.144;
my $global_counter = 0;
my $helix_score = 220;
my $re_helix_score = 178;
my $signal = 1;

my %models = ('MEMSAT-SVM_w27_RE.model'=>27, 'MEMSAT-SVM_w27_SP.model'=>27,
              'MEMSAT-SVM_w33_HL.model'=>33, 'MEMSAT-SVM_w35_IO.model'=>35,
              'MEMSAT-SVM_w33_GM.model'=>33);
my @windows = ("27", "33", "35", "GM");

my (@mtx, $blast_out, $filename, $svm_all, $memsat_out, $memsat_out_single, $globmem_out,
    $png_schematic_out, $png_cartoon_out, $png_cartoon_memsat3_out,
    %range, %raw_hl, %raw_io, %raw_re, %raw_sp, %mtx_list, $header, @constraints);
my ($HL_prediction, $IO_prediction, $RE_prediction, $SP_prediction, $GM_prediction, $system);
my %topology = ();

&main();
exit;

# Subroutines

sub main {

	print "\nMEMSAT-SVM: Alpha-helical transmembrane protein topology prediction\n";
	print "using Support Vector Machines\n\n";

	&get_arguments();

	# If we've been passed .mtx files rather than fasta files
	if ($mtx){

		foreach (@ARGV){
			if ($_ =~ /constraints/){
				push @constraints,$_;
			}else{
				push @mtx,$_;
			}
		}

	# Otherwise run PSI-BLAST to create .mtx files
	}else{

		my @fastas;
		foreach (@ARGV){

			# Deal with constraints
			if ($_ =~ /constraints/){
				push @constraints,$_;
			}else{
				push @fastas,$_;
				&run_psiblast($_);
			}
		}
		@ARGV = @fastas;
	}

	# Create filenames for input files if we are processing multiple sequences
	if ($mtx[0] ne $mtx[-1]){

		my ($x,$y);
		if ($mtx[0] =~ /\//){
			my @tmp = split(/\//,$mtx[0]);
			$x = $tmp[-1];
		}else{
			$x = $mtx[0];
		}
		if ($mtx[-1] =~ /\//){
			my @tmp = split(/\//,$mtx[-1]);
			$y = $tmp[-1];
		}else{
			$y = $mtx[-1];
		}

		$header = $x."-".$y;
	}else{
		if ($mtx[0] =~ /\//){
			my @tmp = split(/\//,$mtx[0]);
			$header = $tmp[-1];
		}else{
			$header = $mtx[0];
		}
	}
	$header =~ s/\.mtx//g;
	$header =~ s/\.fasta//g;
	$header =~ s/\.fa//g;

	# Create filenames for SVM prediction files if we are processing multiple sequences
	$HL_prediction = $output_path.$header."_SVM_w33_HL.prediction";
	$IO_prediction = $output_path.$header."_SVM_w35_IO.prediction";
	$RE_prediction = $output_path.$header."_SVM_w27_RE.prediction";
	$SP_prediction = $output_path.$header."_SVM_w27_SP.prediction";
	$GM_prediction = $output_path.$header."_SVM_w33_GM.prediction";

	# Remove input files from previous runs
	foreach my $req_win (@windows){
		my $input = $input_path.$header."_w".$req_win.".input";
		my $erase = `rm $input` if -e $input;
	}

	print "Generating SVM input files...\n";
	# Create SVM classify input files
	if($runmem3 != 2)
	{
		foreach (@mtx){

			my @tmp = split(/\./,$_);
			$tmp[0] =~ s/output\///;
			&create_input($_);
		}
	}
	# Do SVM classification
	if($runmem3 != 2)
	{
		&classify();
	}
	# Parse predictions and run memsat-svm, globmem-svm and produce graphics
	foreach (@mtx){

		my @tmp2 = split(/\//, $_);
		my @tmp = split(/\./, $tmp2[-1]);

		$filename = (defined($memsat_out_single) && (scalar(@mtx) == 1)) ? $memsat_out_single : $tmp[0];

		$memsat_out = $output_path.$filename.".memsat_svm";
		$globmem_out = $output_path.$filename.".globmem_svm";
		$svm_all = $output_path.$filename."_SVM_ALL.out";
		$png_schematic_out = $output_path.$filename."_schematic.png";
		$png_cartoon_out = $output_path.$filename."_cartoon_memsat_svm.png";
		$png_cartoon_memsat3_out = $output_path.$filename."_cartoon_memsat3.png";

		&parse_predictions($_) if ($runmem3 != 2);
		&run_memsat3($_) if ($runmem3);
		&run_memsat($_) unless (($globmem == 2) || ($runmem3 == 2));

		print "Written file $globmem_out\n" if ($globmem);
		print "\n";
	}

	# Clean up
	$system = `rm $input_path/$header* &> /dev/null` if $remove_files;
	$system = `rm $HL_prediction $IO_prediction $RE_prediction $SP_prediction &> /dev/null` if $remove_files;
}

# Process command line arguments
sub get_arguments {

	my $result = GetOptions ("n=s" => \$ncbidir,
				"mtx=i" => \$mtx,
				"d=s" => \$dbname,
				"o=s" => \$memsat_out_single,
				"i=s" => \$input_path,
				"j=s" => \$output_path,
				"w=s" => \$mem_dir,
				"e=i" => \$remove_files,
				"p=i" => \$globmem,
				"g=i" => \$graphics,
				"m=i" => \$helix_score,
				"r=i" => \$re_helix_score,
				"f=i" => \$format,
				"3=i" => \$runmem3,
				"s=i" => \$signal,
				"c=i" => \$cores,
			        "h"  => sub {&usage;});

	&usage if (!$ARGV[0]);

	## Make paths absolute.
	my $currentWD = cwd();
	foreach my $path ($mem_dir, $input_path, $output_path, $ncbidir, $dbname)
	{
		next unless (defined($path));
		$path = "$currentWD/$path" unless ((substr($path, 0, 1) eq '/') || (substr($path, 0, 1) eq '~'));
	}

	## Make sure directories end with slash...
	foreach my $path ($mem_dir, $input_path, $output_path)
	{
		next unless (defined($path));
		$path .= '/' unless ($path =~ m{/$});
	}

	## ...but get rid of any trailing slashes to the NCBI directory.
	$ncbidir =~ s/\/+$//;

	unless($mtx){

		## Check the NCBI directory
		unless (-d $ncbidir){

			## Look for the NCBI directory
			my $system = `which blastpgp`;
			if ($system =~ /(.*)\/blastpgp$/){
				$ncbidir = $1;
				$ncbidir =~ s/\s+//g;
			}

			unless (-d $ncbidir){
				print "NCBI directory $ncbidir doesn't exist. Please pass it using\n";
				print "the -n paramater or modifiy the value at the top of the script.\n\n";
				exit 1;
			}
		}

		## Make sure we can find blastpgp & makemat
		my $psiblast = $ncbidir."/blastpgp";
		my $makemat = $ncbidir."/makemat";
		unless (-e $psiblast){
			print "Can't find the program blastpgp in the NCBI directory $ncbidir\n";
			print "Please pass the correct NCBI location using the -n parameter or modify\n";
			print "the value at the top of the script..\n\n";
			exit 1;
		}
		unless (-e $makemat){
			print "Can't find the program makemat in the NCBI directory $ncbidir\n";
			print "Please pass the correct NCBI location using the -n parameter or modify\n";
			print "the value at the top of the script.\n\n";
			exit 1;
		}

		unless (-T $dbname){
      # print $dbname
			print "The database name for PSI-BLAST searches has not been set correctly.\n";
			print "Please pass it using the -d parameter or modify the value at the top of the script.\n\n";
			exit 1;
		}
	}

	## Check that $mem_dir is valid and then set all other paths.
	my $this_exe = $0;
	$this_exe =~ s{^.*/}{};

	unless (-e $mem_dir.$this_exe){
		print "The path to the main MEMSAT-SVM directory seems to be incorrect.\n";
		print "Please pass it using the -w parameter or modify the value at the top of the script.\n\n";
		exit 1;
	}

	$input_path = $mem_dir.'input/' unless (defined($input_path));
	$output_path = $mem_dir.'output/' unless (defined($output_path));

	unless ((-d $input_path) && (-w $input_path) && (-d $output_path) && (-w $output_path)){
		print "Cannot find (or write to) the 'input' and 'output' directories.\n";
		print "Please pass valid values using the '-i' and '-j' flags or use the default folders.\n\n";
		exit 1;
	}

	$model_path = $mem_dir.'models/';
	$datadir = $mem_dir.'data/';
	$svm_classify = $mem_dir.'bin/svm_classify';
	$memsat_svm_bin_path = $mem_dir.'bin/memsat-svm';
	$globmem_bin_path = $mem_dir.'bin/globmem';
	$mem_pred_svm_bin_path = $mem_dir.'bin/mem_pred';
	$nnsat_bin_path = $mem_dir.'bin/nnsat';

	unless (-e $memsat_svm_bin_path){
		print "Can't find $memsat_svm_bin_path - have you run make?\n\n";
		exit 1;
	}

	unless (-e $svm_classify){
		print "Can't find $svm_classify - have you run make?\n\n";
		exit 1;
	}

}

# Usage
sub usage {

	print "Version 1.3\n\n";
	print "Usage: run_memsat-svm.pl [options] <fasta file 1> [<constraints file 1>] [<fasta file 2> [<fasta file 3>...]]\n\n";
	print "Options:\n\n";
	print "-p <0|1|2>     Programs to run. memsat-svm predicts topology, globmem-svm\n";
	print "               discriminates between transmembrane and globular proteins. Default 0.\n";
	print "               0 = Run memsat-svm\n";
	print "               1 = Run memsat-svm and globmem-svm\n";
	print "               2 = Run globmem-svm\n";
	print "-s <0|1>       Run memsat-svm with signal peptide function. Default: 1\n";
	print "-3 <0|1|2>     Run memsat version 3. Default 0.\n";
	print "               0 = Run memsat-svm\n";
	print "               1 = Run memsat-svm and memsat3\n";
	print "               2 = Run memsat3\n";
	print "-mtx <0|1>     Process PSI-BLAST .mtx files instead of fasta files. Default 0.\n";
	print "-n <directory> NCBI binary directory (location of blastpgp and makemat)\n";
	print "-d <path>      Database for running PSI-BLAST.\n";
	print "-i <path>      Path to folder for SVM classify input files. Default input/\n";
	print "-j <path>      Output path for all files. Default output/\n";
	print "-w <path>      Main MEMSAT-SVM directory that contains run_memsat-svm.pl. Default ''\n";
	print "-o <filename>  Output filename when running a single sequence. Default:\n";
	print "               <fasta file>.memsat_svm\n";
  	print "-f <1|2|3>     Output format for topology string.       Default: 1\n";
        print "               1 = 6-21,40-57,142-172,214-236,276-302\n";
        print "               2 = A.6,21;B.40,57;C.142,172;D.214,236;E.276,302\n";
        print "               3 = i6-21o40-57i142-172o214-236i276-302o\n";
	print "-e <0|1>       Erase intermediate files. Default 0.\n";
	print "-g <0|1>       Draw topology schematic and cartoon. Default 1.\n";
        print "-m <int>       Minimum score for a transmembrane helix. Default: 220\n";
	print "-r <int>       Minimum score for a re-entrant helix.    Default: 178\n";
	print "-h <0|1>       Show help. Default 0.\n";
	print "-c <int>       Number of CPU cores to use for PSI-BLAST. Default 1.\n\n";
	print "A contraints file for fasta file XYZ.fa must be named XYZ.constraints\n\n";
	print "The constraints file should have the following format, where s,o,m,i\n";
	print "are signal peptide, outside loop, membrane and inside loop:\n";
        print "s:   1-15\n";
        print "o:   1-30\n";
        print "m:   37-59,82-100\n";
        print "i:   65,80,220-230\n\n";
	exit 1;
}

# Run PSI-BLAST
sub run_psiblast {

	my $fasta = shift;
	my $mtx;

	if ($fasta =~ /\.mtx$/){
		print "This looks like an .mtx file! It should be a fasta file, otherwise\n";
		print "pass the -mtx 1 flag.\n\n";
		exit 1;
	}

	if ($fasta =~ /\//){
		my @tmp = split(/\//,$fasta);
		$mtx = $tmp[-1].".mtx";
	}else{
		$mtx = $fasta.".mtx";
	}

	my $out_mtx = $output_path.$mtx;

        my $tmp_rootname = 'memsat-svm_tmp';
	my $tmp_rootpath = $output_path.$tmp_rootname;
	my $blast_out = "$tmp_rootpath.out";

	die "Fasta file $fasta doesn't exist!\n" unless -e $fasta;

	unless (-e $out_mtx || -e $mtx){

		print "Running PSI-BLAST: $fasta\n";

		my $tmp_fasta = "$tmp_rootpath.fasta";
		my $tmp_chk = "$tmp_rootpath.chk";
		my $system = `cp -f $fasta $tmp_fasta`;
		print "$ncbidir/blastpgp -a $cores -j 2 -h 1e-3 -e 1e-3 -b 0 -d $dbname -i $tmp_fasta -C $tmp_chk >& $blast_out\n\n";
		$system = `$ncbidir/blastpgp -a $cores -j 2 -h 1e-3 -e 1e-3 -b 0 -d $dbname -i $tmp_fasta -C $tmp_chk >& $blast_out`;

		unless (-e $tmp_chk){

			print "There was an error running PSI-BLAST. Did you set the database path correctly?\n\n";
			open(ERROR,$blast_out);
			my @error = <ERROR>;
			close ERROR;
			foreach my $line (@error){
				print $line;
			}
			print "\n";
			exit 1;
		}

		$system = `echo $tmp_rootname.chk > $tmp_rootpath.pn`;
		$system = `echo $tmp_rootname.fasta > $tmp_rootpath.sn`;
		$system = `echo "$ncbidir/makemat -P $tmp_rootpath"`;
		$system = `$ncbidir/makemat -P $tmp_rootpath`;
		$system = `cp $tmp_rootpath.mtx $out_mtx`;
		$system = `rm -f ${tmp_rootpath}*`;
		unlink('error.log');

	}

	if (-e $mtx){
		$system = `mv $mtx $output_path`;
		$mtx = $output_path.$mtx;
		push @mtx,$mtx;
	}elsif (-e $out_mtx){
		$mtx = $out_mtx;
		push @mtx,$mtx;
	}else{
		die "Problem creating $mtx\n";
	}
}

# Run memsat-svm and draw images
sub run_memsat {

	my $mtx = shift;
	my $seq = $mtx_list{$mtx};
	my $constrained_prediction = "";

	print "Running MEMSAT-SVM...\n";

	foreach my $c (@constraints){
		my ($x,$y) = split(/\./,$c);
		if ($x =~ /$filename/){
			$constrained_prediction = "-c ".$c." ";
		}
	}

	my $run_command = "$memsat_svm_bin_path ".$constrained_prediction."-m $helix_score -r $re_helix_score -s $signal -f $format $svm_all > $memsat_out";
	print "$run_command\n";
	my $system = `$run_command`;

	if (-e $memsat_out){

		$system = `cat $memsat_out | grep "sequence length is"`;
		if ($system){
			print "Sequence length must be between 30 and 2000 residues.\n"
		}else{

			$system = `cat $memsat_out | grep "No transmembrane helices predicted"`;

			if ($system){
				print "No transmembrane helices predicted.\n"
			}else{

				## check signal score
				$system = `cat $memsat_out | grep "^Signal peptide"`;
				my @split = split(/\s+/,$system);
				my $pred_signal = 'Not detected.';
				if ($split[2] ne 'Not'){
					if ($format == 1 || $format == 3){

						my ($start,$stop) = split(/-/,$split[2]);
						if ($format == 3){
							$stop =~ s/o//g;
						}
						$pred_signal = $stop;
					}else{
						my ($start,$stop) = split(/,/,$split[2]);
						$pred_signal = $stop;
					}
				}

				## check predicted re-entrant region
				$system = `cat $memsat_out | grep "^Re-entrant helices:"`;
				@split = split(/\s+/,$system);
				my @re_helix_split = split(/\D+/,$split[2]);
				shift @re_helix_split if $format == 3;
				shift @re_helix_split if $format == 2;
				my $re_helix_array = \@re_helix_split;

				$system = `cat $memsat_out | grep "^N-terminal:"`;
				@split = split(/\s+/,$system);
				my $pred_n = $split[1];

				$system = `cat $memsat_out | grep "^Topology:"`;
				@split = split(/\s+/,$system);
				my @topology_split = split(/\D+/,$split[1]);
				shift @topology_split if $format == 3;
				shift @topology_split if $format == 2;
				my $topology_array = \@topology_split;

				print "\nWritten file $memsat_out\n";

				unless(`perl -MGD -e 1 &> /dev/stdout`){

					#print "@$topology_array\n$pred_n\n@$re_helix_array\n$pred_signal\n";
					&draw_image($filename,$seq,$topology_array,$pred_n, $re_helix_array,$pred_signal) if $graphics;
					print "Written file $png_schematic_out\n" if $graphics;
					print "Written file $png_cartoon_out\n" if $graphics;
					print "Written file $png_cartoon_memsat3_out\n" if $graphics && -e $png_cartoon_memsat3_out && $runmem3;

				}else{
					print "Couldn't find GD module. Graphics disabled.\n" if $graphics;
				}
			}
		}

	}else{

		die "Problem running $memsat_svm_bin_path\n";
	}

}

sub run_memsat3 {

	my $mtx = shift;
	my $seq = $mtx_list{$mtx};

	print "\nRunning MEMSAT3...\n";

	my $rootname = $output_path.$filename;

	my $globmem = "$rootname.globmem";
	my $nn = "$rootname.nn";
	my $memsat3 = "$rootname.memsat3";

	my $run_command ="$globmem_bin_path ".$datadir."glob_weights.dat $mtx > $globmem";
	print "$run_command\n";
	my $system = `$run_command`;
	$run_command = "$mem_pred_svm_bin_path ".$datadir."weights.dat $mtx > $nn";
	print "$run_command\n";
	$system = `$run_command`;
	$run_command = "$nnsat_bin_path $nn > $memsat3";
	print "$run_command\n";
	$system = `$run_command`;

	print "\nWritten file $globmem\n";
	print "Written file $nn\n";
	print "Written file $memsat3\n";
	print "\n" unless $runmem3 == 2;

	unless (open (MEM3,$memsat3)){
			die "Can't open $memsat3!\n";
	}
	my @memsat3 = <MEM3>;
	close MEM3;

	my $tag = 0;
	foreach my $line (@memsat3){

		if ($line =~ /^================/){
			$tag++;
			next;
		}
		next unless $tag;

		if ($line =~ /^\d+:\s+\((\w+)\)\s(\d+)-(\d+)\s+\(-*\d+\.\d+\)/){
			$topology{'MEMSAT3'}{'n_term'} = $1;
			push @{$topology{'MEMSAT3'}{'topology_array'}},$2;
			push @{$topology{'MEMSAT3'}{'topology_array'}},$3;
		}

		if ($line =~ /^\d+:\s+(\d+)-(\d+)\s+\(-*\d+\.\d+\)/){
			push @{$topology{'MEMSAT3'}{'topology_array'}},$1;
			push @{$topology{'MEMSAT3'}{'topology_array'}},$2;
		}
	}

	if ($graphics && $runmem3 == 2){
		&draw_image($filename,$seq,\%topology,'',[],0);
		print "Written file $png_schematic_out\n";
		print "Written file $png_cartoon_memsat3_out\n";
	}
}

# Parse SVM output
sub parse_predictions {

	my $mtx = shift;

	print "Parsing SVM output files...\n";

	my $seq = $mtx_list{$mtx};
	my @seq = split(//,$seq);
	my $length = length $seq;

	if ($globmem != 2){


		unless (open (HL,$HL_prediction)){
			die "Can't open $HL_prediction!\n";
		}
		my @hl = <HL>;
		close HL;

		for($global_counter..($global_counter+$length-1)){
			$hl[$_] =~ s/\s+//g;
			$hl[$_] = Round::nearest_ceil(.001,$hl[$_]);
		}

			unless (open (IO,$IO_prediction)){
		die "Can't open $IO_prediction!\n";
		}
		my @io = <IO>;
		close IO;

		for($global_counter..($global_counter+$length-1)){
			$io[$_] =~ s/\s+//g;
			$io[$_] = Round::nearest_ceil(.001,$io[$_]);
		}

		unless (open (RE,$RE_prediction)){
			die "Can't open $RE_prediction!\n";
		}
		my @re = <RE>;
		close RE;

		for($global_counter..($global_counter+$length-1)){
			$re[$_] =~ s/\s+//g;
			$re[$_] = Round::nearest_ceil(.001,$re[$_]);
		}

		unless (open (SP,$SP_prediction)){
			die "Can't open $SP_prediction!\n";
		}
		my @sp = <SP>;
		close SP;

		for($global_counter..($global_counter+$length-1)){
			$sp[$_] =~ s/\s+//g;
			$sp[$_] = Round::nearest_ceil(.001,$sp[$_]);
		}

		if (scalar @hl != scalar @io || scalar @hl != scalar @re || scalar @re != scalar @sp){
			die "Uneven number of values for $header\n";
		}

		open (RAW,">$svm_all");
		my $seq_counter = 0;
		for ($global_counter..($global_counter+$length-1)){

			print RAW "$seq[$seq_counter]\t$hl[$_]\t$io[$_]\t$re[$_]\t$sp[$_]\n";

			#print "$seq[$_]\t$hl[$_]\t$io[$_]\t$re[$_]\t$sp[$_]\n";

			$raw_hl{$seq_counter+1} = $hl[$_];
			$raw_io{$seq_counter+1} = $io[$_];
			$raw_re{$seq_counter+1} = $re[$_];
			$raw_sp{$seq_counter+1} = $sp[$_];

			$seq_counter++;

		}
		close RAW;

	}

	if ($globmem){

		print "Running GLOBMEM-SVM...\n";

		unless (open (GM,$GM_prediction)){
			die "Can't open $GM_prediction!\n";
		}
		my @gm = <GM>;
		close GM;

		my $tm_residues = 0;
		my $tm_score = 0;

		for($global_counter..($global_counter+$length-1)){
			$gm[$_] =~ s/\s+//g;
			$tm_residues++ if $gm[$_] > $globmem_score_threshold;
			if ($gm[$_] > 0){
				$tm_score += $gm[$_];
			}
		}

		if (1){

			open (GMO,">$globmem_out");
			print GMO "Transmembrane residues found:\t$tm_residues\n";
			print GMO "Transmembrane score:\t\t$tm_score\n";
			# Using the same threshold as in the benchmarking, see the publication.
			if ($tm_residues > 0){
				print GMO "This looks like a transmembrane protein.\n";
			}else{
				print GMO "This looks like a globular protein.\n";
			}
			close GMO;

		}
	}

	$global_counter += $length;

}

# Run SVM classify
sub classify {

	print "Running svm_classify...\n";

	foreach my $model (keys %models){

		my $type;
		if ($model =~ /HL/){
			$type = 'HL';
		}elsif($model =~ /IO/){
			$type = 'IO';
		}elsif($model =~ /SP/){
			$type = 'SP';
		}elsif($model =~ /RE/){
			$type = 'RE';
		}elsif($model =~ /GM/){
			$type = 'GM';
			next unless $globmem;
		}

		if ($globmem == 2){
			next unless $type eq 'GM';
		}

		my $input;

		if ($type eq 'GM'){
			$input = $input_path.$header."_w".$models{$model}."_GM.input";
		}else{
			$input = $input_path.$header."_w".$models{$model}."_TM.input";
		}
		my $prediction = $output_path.$header."_SVM_w".$models{$model}."_".$type.".prediction";;
		$model = $model_path.$model;

		if (-e $model){

			print "$svm_classify -v 0 $input $model $prediction\n" unless -e $prediction;
			my $system = `$svm_classify -v 0 $input $model $prediction` unless -e $prediction;

		}else{
			die "$model doesn't exist.\n";
		}
	}
	print "\n";
}

# Create input files for SVM classify
sub create_input {

	my $mtx = shift;
	my @mtx = ();

	if (-e $mtx){
		open (MTX,$mtx);
		@mtx = <MTX>;
		close MTX;
		$mtx[1] =~ s/\s+//g;
		$mtx_list{$mtx} = $mtx[1];

	}else{
		die "Couldn't find $mtx\n";
	}

	if ($mtx[0] =~ /\>/){
		print "This looks like a fasta file! It should be a .mtx file, otherwise\n";
		print "pass the -mtx 0 flag.\n\n";
		exit 1;
	}

	my $length = length $mtx_list{$mtx};
	my $seq = $mtx_list{$mtx};
	my @empty = ();

	for (1..20){
		push @empty,0;
	}

	foreach my $reqwin (@windows){

		if ($globmem == 2){
			next unless $reqwin eq "GM";
		}

		my ($input,$win);

		next if ($globmem == 0 && $reqwin eq "GM");

		if ($reqwin eq "GM"){
			&get_gm_normalisation_values();
			$win = 33;
			$input = $input_path.$header."_w".$win."_GM.input";
			#print "Creating $mtx $win $input\n";
		}else{
			get_normalisation_values();
			$win = $reqwin;
			$input = $input_path.$header."_w".$win."_TM.input";
			#print "Creating $mtx $win $input\n";
		}


		## Add X's to either side of seqeunce
		my $x_string;
		for (1..($win - 1)/2){
			$x_string .= 'X';
		}

		my $seq_x = $x_string.$seq.$x_string;

		## Fill up %profile with lines from @mtx

		my %profile = ();

		for (1..($length + $win - 1)){

			if (($_ <= ($win - 1)/2)||($_ > $length + (($win - 1)/2))){

				## Missing
				@{$profile{$_}} = @empty;

			}else{
				## Start on MTX line 15

				if (defined $mtx[$_+13-(($win - 1)/2)]){

					@{$profile{$_}} = split(/\s+/,$mtx[$_+13-(($win - 1)/2)]);

					## Now do the Z score normalisation

					## Not normalised
					#print "Not normalised: @{$profile{$_}}\n";

					my $z_count = 0;
					foreach my $z (@{$profile{$_}}){

						## Normalise to Z score
						$z= ($z - $range{$z_count}{'mean'})/$range{$z_count}{'sd'};
						## Scale
						$z = ($z + (-1 * $range{$z_count}{'lower'}))/$range{$z_count}{'range'};
						## 5 Decimal places
						$z = Round::nearest_ceil(.00001,$z);
						$z_count++;
					}

    					## Get 20 residues rather than 28
					my @aa_20 = (${$profile{$_}}[1],${$profile{$_}}[3],${$profile{$_}}[4],${$profile{$_}}[5],${$profile{$_}}[6],
					${$profile{$_}}[7],${$profile{$_}}[8],${$profile{$_}}[9],${$profile{$_}}[10],${$profile{$_}}[11],
					${$profile{$_}}[12],${$profile{$_}}[13],${$profile{$_}}[14],${$profile{$_}}[15],${$profile{$_}}[16],
					${$profile{$_}}[17],${$profile{$_}}[18],${$profile{$_}}[19],${$profile{$_}}[21],${$profile{$_}}[22]);

					@{$profile{$_}} = @aa_20;

				}
			}
		}

		open (IN,">>$input");

		for my $pos ((($win - 1)/2) + 1..($length + ($win - 1)/2)){

			my $real_pos = $pos - (($win - 1)/2);
			my ($vector,$seq,@array);

			## Push all the arrays onto @array

			for my $w ($pos - (($win - 1)/2) .. $pos + (($win - 1)/2)){
				push @array,@{$profile{$w}}
			}

			## Generate sequence window for the comment

			for my $j ($real_pos - 1 - (($win - 1)/2) .. $real_pos - 1 + (($win - 1)/2)){
				if ($j < 0 || $j >= $length){
					$seq .= 'X';
				}else{
					$seq .= substr($seq_x,$j+(($win - 1)/2),1);
				}
			}

			my $count = 1;
			foreach my $v (@array){
				$vector .= $count.":".$v." ";
				$count++;
			}
			print IN "0 $vector # $seq\n";
		}

		close IN;
	}
}

# Draw images
sub draw_image{

	my ($header,$seq,$topology,$in_out,$re_region,$signal_region) = @_;

	load_module('DrawTransmembraneSchematic');
	load_module('DrawTransmembraneCartoon');

 	my $svm_raw_muliplier = 6;
	my $svm_raw_sp_muliplier = 6;
	my $svm_raw_prob = 40;

	$topology{'MEMSAT-SVM'}{'topology_array'} = $topology;
	$topology{'MEMSAT-SVM'}{'n_term'} = $in_out;
	if (scalar @$re_region){
		$topology{'MEMSAT-SVM'}{'reentrant_array'} = $re_region;
	}
	if ($signal_region !~ 'detected'){
		$topology{'MEMSAT-SVM'}{'signal'} = $signal_region;
	}

	my @order = ('MEMSAT-SVM');
	push @order,'MEMSAT3' if $runmem3 == 1;

	if ($runmem3 == 2){

		@order = ('MEMSAT3') ;

		my $im = DrawTransmembraneSchematic->new(-title=>$header,
							 -inside_rgb=>,[256,5,50],
						 	 -outside_rgb=>,[140,399,70],
						 	 -signal_rgb=>,[633,23,46],
						 	 -sequence=>$seq,
                                                 	 -topologies=>\%topology,
						 	 -draw_consensus=>0,
						 	 -draw_key=>1,
						 	 -draw_custom_plot_1=>0,
						 	 -draw_custom_plot_3=>0,
						 	 -draw_custom_plot_5=>0,
						 	 -draw_custom_plot_7=>0,
						 	 -order=>\@order
					         	);

		open(OUTPUT, ">$png_schematic_out");
		binmode OUTPUT;
		print OUTPUT $im->png;
		close OUTPUT;

	}else{

		my $im = DrawTransmembraneSchematic->new(-title=>$header,
							 -inside_rgb=>,[256,5,50],
						 	 -outside_rgb=>,[140,399,70],
						 	 -signal_rgb=>,[633,23,46],
						 	 -sequence=>$seq,
                                                 	 -topologies=>\%topology,
						 	 -draw_consensus=>0,
						 	 -draw_key=>1,
						 	 -unknown_label=>'Re-entrant Helix',
						 	 -draw_custom_plot_1=>1,
						 	 -custom_plot_multiplier_1=>$svm_raw_muliplier,
						 	 -custom_plot_label_1=>'SVM H/L Raw',
						 	 -custom_plot_data_1=>\%raw_hl,
						 	 -draw_custom_plot_3=>1,
						 	 -custom_plot_multiplier_3=>$svm_raw_muliplier,
						 	 -custom_plot_label_3=>'SVM iL/oL Raw',
						 	 -custom_plot_data_3=>\%raw_io,
						 	 -draw_custom_plot_5=>1,
						 	 -custom_plot_multiplier_5=>$svm_raw_muliplier,
						 	 -custom_plot_label_5=>'SVM RE/!RE Raw',
						 	 -custom_plot_data_5=>\%raw_re,
						 	 -draw_custom_plot_7=>1,
						 	 -custom_plot_multiplier_7=>$svm_raw_sp_muliplier,
							 -custom_plot_label_7=>'SVM SP/!SP Raw',
						 	 -custom_plot_data_7=>\%raw_sp,
						 	 -order=>\@order
					         	);

		open(OUTPUT, ">$png_schematic_out");
		binmode OUTPUT;
		print OUTPUT $im->png;
		close OUTPUT;

	}


	my %labels;

	if (scalar @$re_region){
		for (my $r = 0; $r < scalar @$re_region; $r += 2){
			$labels{$$re_region[$r]} = $$re_region[$r]."-".$$re_region[$r+1]." Re-entrant helix";
		}
	}

	if ($signal_region !~ 'detected'){
		$labels{1} = "1-".$signal_region." Signal peptide";
	}

	if ($runmem3 != 2){

	 	my $im = DrawTransmembraneCartoon->new(-title=>$header,
	                                     	       -n_terminal=>$in_out,
					     	       -topology=>$topology,
					     	       -labels=> \%labels,
					     	       -loop_width=>25,
					     	       -text_offset=>-25,
					     	       -outside_label=>'Extracellular',
					     	       -inside_label=>'Cytoplasmic',
					     	       -membrane_label=>'Membrane',
					     	       -colour_scheme=>'yellow');

		open(OUTPUT, ">$png_cartoon_out");
		binmode OUTPUT;
		print OUTPUT $im->png;
		close OUTPUT;

	}

	if ($runmem3){

		my $im = DrawTransmembraneCartoon->new(-title=>$header,
	                                     	       -n_terminal=>$topology{'MEMSAT3'}{'n_term'},
					               -topology=>$topology{'MEMSAT3'}{'topology_array'},
					               -loop_width=>25,
					               -text_offset=>-25,
					               -outside_label=>'Extracellular',
					               -inside_label=>'Cytoplasmic',
					               -membrane_label=>'Membrane',
					               -colour_scheme=>'yellow');

		open(OUTPUT, ">$png_cartoon_memsat3_out");
		binmode OUTPUT;
		print OUTPUT $im->png;
		close OUTPUT;
	}

	%raw_hl = ();
	%raw_io = ();
	%raw_re = ();
	%raw_sp = ();
	%topology = ();

}

# Load graphics modules
sub load_module {
	eval "require $_[0]";
    die if $@;
    $_[0]->import(@_[1 .. $#_]);
}

# Normalisation values
sub get_normalisation_values {

        $range{0}{'mean'} = -32768;
        $range{0}{'sd'} = 1;
        $range{0}{'lower'} = 0;
        $range{0}{'upper'} = 0;
        $range{0}{'range'} = 1;
        $range{1}{'mean'} = -65.7855796773707;
        $range{1}{'sd'} = 165.265945369071;
        $range{1}{'lower'} = -2.67577460882869;
        $range{1}{'upper'} = 4.54289344365968;
        $range{1}{'range'} = 7.21866805248836;
        $range{2}{'mean'} = -32768;
        $range{2}{'sd'} = 1;
        $range{2}{'lower'} = 0;
        $range{2}{'upper'} = 33266;
        $range{2}{'range'} = 33266;
        $range{3}{'mean'} = -263.751620684456;
        $range{3}{'sd'} = 123.206910338565;
        $range{3}{'lower'} = -3.37844994385232;
        $range{3}{'upper'} = 11.223003781888;
        $range{3}{'range'} = 14.6014537257403;
        $range{4}{'mean'} = -242.89740690487;
        $range{4}{'sd'} = 188.205050583637;
        $range{4}{'lower'} = -2.4234343960522;
        $range{4}{'upper'} = 5.73787687182693;
        $range{4}{'range'} = 8.16131126787913;
        $range{5}{'mean'} = -194.765716870194;
        $range{5}{'sd'} = 180.643157809653;
        $range{5}{'lower'} = -2.38168048182125;
        $range{5}{'upper'} = 5.31858349089852;
        $range{5}{'range'} = 7.70026397271977;
        $range{6}{'mean'} = -135.345507311925;
        $range{6}{'sd'} = 251.687254013587;
        $range{6}{'lower'} = -2.19182530656994;
        $range{6}{'upper'} = 4.06196774373301;
        $range{6}{'range'} = 6.25379305030295;
        $range{7}{'mean'} = -192.372795115332;
        $range{7}{'sd'} = 230.165335412046;
        $range{7}{'lower'} = -2.14466354805761;
        $range{7}{'upper'} = 4.12908743801041;
        $range{7}{'range'} = 6.27375098606802;
        $range{8}{'mean'} = -221.253316749585;
        $range{8}{'sd'} = 155.427737434726;
        $range{8}{'lower'} = -2.65555357147081;
        $range{8}{'upper'} = 8.40424841993316;
        $range{8}{'range'} = 11.059801991404;
        $range{9}{'mean'} = -92.5749283883612;
        $range{9}{'sd'} = 246.71482166236;
        $range{9}{'lower'} = -2.51474584068856;
        $range{9}{'upper'} = 3.33411232793338;
        $range{9}{'range'} = 5.84885816862194;
        $range{10}{'mean'} = -179.119704507764;
        $range{10}{'sd'} = 167.0612436708;
        $range{10}{'lower'} = -2.51931738471669;
        $range{10}{'upper'} = 5.86084290403797;
        $range{10}{'range'} = 8.38016028875466;
        $range{11}{'mean'} = -89.1816674204734;
        $range{11}{'sd'} = 236.538032097989;
        $range{11}{'lower'} = -2.59078139419735;
        $range{11}{'upper'} = 3.07426954122637;
        $range{11}{'range'} = 5.66505093542372;
        $range{12}{'mean'} = -88.3949193426805;
        $range{12}{'sd'} = 185.952348451384;
        $range{12}{'lower'} = -2.88033512412105;
        $range{12}{'upper'} = 5.65410938930711;
        $range{12}{'range'} = 8.53444451342817;
        $range{13}{'mean'} = -188.69587667722;
        $range{13}{'sd'} = 183.448068516649;
        $range{13}{'lower'} = -2.52008171610062;
        $range{13}{'upper'} = 5.73838626478471;
        $range{13}{'range'} = 8.25846798088533;
        $range{14}{'mean'} = -249.184494195688;
        $range{14}{'sd'} = 178.122836782516;
        $range{14}{'lower'} = -2.32320298328491;
        $range{14}{'upper'} = 6.35058656502781;
        $range{14}{'range'} = 8.67378954831272;
        $range{15}{'mean'} = -162.548130559325;
        $range{15}{'sd'} = 164.437104518497;
        $range{15}{'lower'} = -2.5690787409429;
        $range{15}{'upper'} = 6.32185864378553;
        $range{15}{'range'} = 8.89093738472843;
        $range{16}{'mean'} = -207.047225991256;
        $range{16}{'sd'} = 173.171717333037;
        $range{16}{'lower'} = -2.25182714599291;
        $range{16}{'upper'} = 5.99432310297486;
        $range{16}{'range'} = 8.24615024896777;
        $range{17}{'mean'} = -90.7070707070707;
        $range{17}{'sd'} = 163.603610103856;
        $range{17}{'lower'} = -2.83791371717405;
        $range{17}{'upper'} = 4.99198685278779;
        $range{17}{'range'} = 7.82990056996184;
        $range{18}{'mean'} = -86.9080732700136;
        $range{18}{'sd'} = 145.736036521422;
        $range{18}{'lower'} = -2.80021288125261;
        $range{18}{'upper'} = 5.76321473616139;
        $range{18}{'range'} = 8.56342761741401;
        $range{19}{'mean'} = -80.3608095884215;
        $range{19}{'sd'} = 211.306064426292;
        $range{19}{'lower'} = -2.69106895703854;
        $range{19}{'upper'} = 3.6599082552985;
        $range{19}{'range'} = 6.35097721233704;
        $range{20}{'mean'} = -278.93061209106;
        $range{20}{'sd'} = 192.05331159166;
        $range{20}{'lower'} = -2.45801222585864;
        $range{20}{'upper'} = 8.10676264412129;
        $range{20}{'range'} = 10.5647748699799;
        $range{21}{'mean'} = -100;
        $range{21}{'sd'} = 1;
        $range{21}{'lower'} = 0;
        $range{21}{'upper'} = 0;
        $range{21}{'range'} = 1;
        $range{22}{'mean'} = -167.07315694256;
        $range{22}{'sd'} = 185.875809713378;
        $range{22}{'lower'} = -2.53355637715103;
        $range{22}{'upper'} = 6.16042054481576;
        $range{22}{'range'} = 8.69397692196679;
        $range{23}{'mean'} = -32768;
        $range{23}{'sd'} = 1;
        $range{23}{'lower'} = 0;
        $range{23}{'upper'} = 33219;
        $range{23}{'range'} = 33219;
        $range{24}{'mean'} = -32768;
        $range{24}{'sd'} = 1;
        $range{24}{'lower'} = 0;
        $range{24}{'upper'} = 32687;
        $range{24}{'range'} = 32687;
        $range{25}{'mean'} = -401.752826775215;
        $range{25}{'sd'} = 3.56447071169638;
        $range{25}{'lower'} = -11.8522991607579;
        $range{25}{'upper'} = 3.29721513397469;
        $range{25}{'range'} = 15.1495142947326;
        $range{26}{'mean'} = -32768;
        $range{26}{'sd'} = 1;
        $range{26}{'lower'} = 0;
        $range{26}{'upper'} = 32687;
        $range{26}{'range'} = 32687;
        $range{27}{'mean'} = -32768;
        $range{27}{'sd'} = 1;
        $range{27}{'lower'} = 0;
        $range{27}{'upper'} = 33113;
        $range{27}{'range'} = 33113;

}

sub get_gm_normalisation_values {

        $range{0}{'mean'} = -32768;
        $range{0}{'sd'} = 1;
        $range{0}{'lower'} = 0;
        $range{0}{'upper'} = 0;
        $range{0}{'range'} = 1;
        $range{1}{'mean'} = -59.4884990869504;
        $range{1}{'sd'} = 194.286939229385;
        $range{1}{'lower'} = -2.68423344863818;
        $range{1}{'upper'} = 3.90396030785909;
        $range{1}{'range'} = 6.58819375649728;
        $range{2}{'mean'} = -31114.7648487615;
        $range{2}{'sd'} = 7151.78256236352;
        $range{2}{'lower'} = -0.23116406809383;
        $range{2}{'upper'} = 4.42543723509149;
        $range{2}{'range'} = 4.65660130318532;
        $range{3}{'mean'} = -253.278286978508;
        $range{3}{'sd'} = 191.082828623659;
        $range{3}{'lower'} = -2.33784331244812;
        $range{3}{'upper'} = 7.40138869183244;
        $range{3}{'range'} = 9.73923200428055;
        $range{4}{'mean'} = -170.379020929906;
        $range{4}{'sd'} = 264.940513149384;
        $range{4}{'lower'} = -2.248130993595;
        $range{4}{'upper'} = 3.8475769855369;
        $range{4}{'range'} = 6.0957079791319;
        $range{5}{'mean'} = -124.718177880789;
        $range{5}{'sd'} = 244.943476773369;
        $range{5}{'lower'} = -2.45069527887291;
        $range{5}{'upper'} = 3.73848771130184;
        $range{5}{'range'} = 6.18918299017475;
        $range{6}{'mean'} = -179.413201760547;
        $range{6}{'sd'} = 268.111735749738;
        $range{6}{'lower'} = -1.94913809639145;
        $range{6}{'upper'} = 4.10803801139325;
        $range{6}{'range'} = 6.0571761077847;
        $range{7}{'mean'} = -180.469957156904;
        $range{7}{'sd'} = 269.268035764037;
        $range{7}{'lower'} = -2.03711532743448;
        $range{7}{'upper'} = 3.51497330335304;
        $range{7}{'range'} = 5.55208863078752;
        $range{8}{'mean'} = -173.214039893243;
        $range{8}{'sd'} = 221.095383730087;
        $range{8}{'lower'} = -2.2107470172397;
        $range{8}{'upper'} = 5.74509525465225;
        $range{8}{'range'} = 7.95584227189195;
        $range{9}{'mean'} = -128.956138502599;
        $range{9}{'sd'} = 270.230396073888;
        $range{9}{'lower'} = -2.27969862179738;
        $range{9}{'upper'} = 3.47094979739488;
        $range{9}{'range'} = 5.75064841919226;
        $range{10}{'mean'} = -122.300972748982;
        $range{10}{'sd'} = 219.14347443512;
        $range{10}{'lower'} = -2.4673288978591;
        $range{10}{'upper'} = 4.1949730655619;
        $range{10}{'range'} = 6.66230196342099;
        $range{11}{'mean'} = -126.702480451374;
        $range{11}{'sd'} = 261.118468960118;
        $range{11}{'lower'} = -2.32575475019877;
        $range{11}{'upper'} = 3.02430725638174;
        $range{11}{'range'} = 5.35006200658051;
        $range{12}{'mean'} = -112.731481481481;
        $range{12}{'sd'} = 210.706829003488;
        $range{12}{'lower'} = -2.54983913459029;
        $range{12}{'upper'} = 5.50875112577517;
        $range{12}{'range'} = 8.05859026036545;
        $range{13}{'mean'} = -146.358471461348;
        $range{13}{'sd'} = 227.975488811242;
        $range{13}{'lower'} = -2.48115063372883;
        $range{13}{'upper'} = 4.5064426742472;
        $range{13}{'range'} = 6.98759330797603;
        $range{14}{'mean'} = -207.458596713021;
        $range{14}{'sd'} = 252.428214816201;
        $range{14}{'lower'} = -1.99082897152718;
        $range{14}{'upper'} = 4.39910648467555;
        $range{14}{'range'} = 6.38993545620273;
        $range{15}{'mean'} = -116.26493655476;
        $range{15}{'sd'} = 206.743223467249;
        $range{15}{'lower'} = -2.5719588508277;
        $range{15}{'upper'} = 4.86722089207517;
        $range{15}{'range'} = 7.43917974290288;
        $range{16}{'mean'} = -149.347444631737;
        $range{16}{'sd'} = 234.557445472634;
        $range{16}{'lower'} = -2.35188678089788;
        $range{16}{'upper'} = 4.27761925275911;
        $range{16}{'range'} = 6.62950603365699;
        $range{17}{'mean'} = -81.3822048976916;
        $range{17}{'sd'} = 184.547968097115;
        $range{17}{'lower'} = -2.87522967916439;
        $range{17}{'upper'} = 4.38033652298082;
        $range{17}{'range'} = 7.25556620214521;
        $range{18}{'mean'} = -84.4269794446786;
        $range{18}{'sd'} = 174.495864546024;
        $range{18}{'lower'} = -2.7655269757298;
        $range{18}{'upper'} = 4.95958446749562;
        $range{18}{'range'} = 7.72511144322541;
        $range{19}{'mean'} = -106.598539120663;
        $range{19}{'sd'} = 239.930284037695;
        $range{19}{'lower'} = -2.39403484717708;
        $range{19}{'upper'} = 3.55769402993144;
        $range{19}{'range'} = 5.95172887710851;
        $range{20}{'mean'} = -285.688778854708;
        $range{20}{'sd'} = 235.662424950478;
        $range{20}{'lower'} = -2.08480932523942;
        $range{20}{'upper'} = 6.75834842652322;
        $range{20}{'range'} = 8.84315775176265;
        $range{21}{'mean'} = -100;
        $range{21}{'sd'} = 1;
        $range{21}{'lower'} = 0;
        $range{21}{'upper'} = 0;
        $range{21}{'range'} = 1;
        $range{22}{'mean'} = -169.542854801704;
        $range{22}{'sd'} = 234.932153737865;
        $range{22}{'lower'} = -2.13021988363798;
        $range{22}{'upper'} = 4.94841951731671;
        $range{22}{'range'} = 7.07863940095469;
        $range{23}{'mean'} = -31113.0723123566;
        $range{23}{'sd'} = 7159.09167487829;
        $range{23}{'lower'} = -0.231164477673982;
        $range{23}{'upper'} = 4.41355883501712;
        $range{23}{'range'} = 4.6447233126911;
        $range{24}{'mean'} = -31177.2046635763;
        $range{24}{'sd'} = 7032.50018195697;
        $range{24}{'lower'} = -0.226206227552628;
        $range{24}{'upper'} = 4.42178512036996;
        $range{24}{'range'} = 4.64799134792259;
        $range{25}{'mean'} = -400.846069204476;
        $range{25}{'sd'} = 6.86956382022376;
        $range{25}{'lower'} = -10.9401313914972;
        $range{25}{'upper'} = 2.01556744603114;
        $range{25}{'range'} = 12.9556988375284;
        $range{26}{'mean'} = -32768;
        $range{26}{'sd'} = 1;
        $range{26}{'lower'} = 0;
        $range{26}{'upper'} = 32687;
        $range{26}{'range'} = 32687;
        $range{27}{'mean'} = -32768;
        $range{27}{'sd'} = 1;
        $range{27}{'lower'} = 0;
        $range{27}{'upper'} = 33113;
        $range{27}{'range'} = 33113;

}
