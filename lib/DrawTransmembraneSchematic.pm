package DrawTransmembraneSchematic;

use strict;
use warnings;
use GD;
use Data::Dumper;
our @ISA = "GD";
our @EXPORT = qw(new add_topology png);
our $VERSION = '1.0';

sub new {

	my $class = shift;
 	my %options = @_;

	my $self = {

		## general paramaters
		'_reentrant_array' => \@{$options{-reentrant}},
		'_topology_array' => \@{$options{-topology}},
		'_topology_string' => exists $options{-topology_string} ? $options{-topology_string} : 0,
		'_n_term' => exists $options{-n_term} ? $options{-n_term} : 'out',
		'_signal' => exists $options{-signal} ? $options{-signal} : 0,
		'_method' => exists $options{-method} ? $options{-method} : 0,
		'_sequence' => exists $options{-sequence} ? $options{-sequence} : 0,
		'_seq_length' => exists $options{-seq_length} ? $options{-seq_length} : 0,
		'_draw_kd' => exists $options{-draw_kd} ? $options{-draw_kd} : 1,
		'_kd_window' => exists $options{-kd_window} ? $options{-kd_window} : 19,
		'_kd_multiplier' => exists $options{-kd_multiplier} ? $options{-kd_multiplier} : 10,
		'_draw_key' => exists $options{-draw_key} ? $options{-draw_key} : 1,
		'_draw_grid' => exists $options{-draw_grid} ? $options{-draw_grid} : 1,
		'_draw_consensus' => exists $options{-draw_consensus} ? $options{-draw_consensus} : 1,
		'_title' => exists $options{-title} ? $options{-title} : '',
		'_inside_label' => exists $options{-inside_label} ? $options{-inside_label} : "Cytoplasmic",
		'_outside_label' => exists $options{-outside_label} ? $options{-outside_label} : "Extracellular",
		'_helix_label' => exists $options{-helix_label} ? $options{-helix_label} : "Transmembrane Helix",
		'_unknown_label' => exists $options{-unknown_label} ? $options{-unknown_label} : "Unknown",
		'_signal_label' => exists $options{-signal_label} ? $options{-signal_label} : "Signal Peptide",
		'_kd_label' => exists $options{-kd_label} ? $options{-kd_label} : "Kyte-Doolittle",
		'_custom_plot_label' => exists $options{-custom_plot_label} ? $options{-custom_plot_label} : "Custom Plot",
		'_consensus_label' => exists $options{-consensus_label} ? $options{-consensus_label} : "Consensus",

		'_draw_custom_plot_1' => exists $options{-draw_custom_plot_1} ? $options{-draw_custom_plot_1} : 0,
		'_custom_plot_multiplier_1' => exists $options{-custom_plot_multiplier_1} ? $options{-custom_plot_multiplier_1} : 10,
		'_custom_plot_data_1' => \%{$options{-custom_plot_data_1}},
		'_custom_plot_label_1' => exists $options{-custom_plot_label_1} ? $options{-custom_plot_label_1} : "Custom Plot 1",

		'_draw_custom_plot_2' => exists $options{-draw_custom_plot_2} ? $options{-draw_custom_plot_2} : 0,
		'_custom_plot_multiplier_2' => exists $options{-custom_plot_multiplier_2} ? $options{-custom_plot_multiplier_2} : 10,
		'_custom_plot_data_2' => \%{$options{-custom_plot_data_2}},
		'_custom_plot_label_2' => exists $options{-custom_plot_label_2} ? $options{-custom_plot_label_2} : "Custom Plot 2",

		'_draw_custom_plot_3' => exists $options{-draw_custom_plot_3} ? $options{-draw_custom_plot_3} : 0,
		'_custom_plot_multiplier_3' => exists $options{-custom_plot_multiplier_3} ? $options{-custom_plot_multiplier_3} : 10,
		'_custom_plot_data_3' => \%{$options{-custom_plot_data_3}},
		'_custom_plot_label_3' => exists $options{-custom_plot_label_3} ? $options{-custom_plot_label_3} : "Custom Plot 3",

		'_draw_custom_plot_4' => exists $options{-draw_custom_plot_4} ? $options{-draw_custom_plot_4} : 0,
		'_custom_plot_multiplier_4' => exists $options{-custom_plot_multiplier_4} ? $options{-custom_plot_multiplier_4} : 10,
		'_custom_plot_data_4' => \%{$options{-custom_plot_data_4}},
		'_custom_plot_label_4' => exists $options{-custom_plot_label_4} ? $options{-custom_plot_label_4} : "Custom Plot 4",

		'_draw_custom_plot_5' => exists $options{-draw_custom_plot_5} ? $options{-draw_custom_plot_5} : 0,
		'_custom_plot_multiplier_5' => exists $options{-custom_plot_multiplier_5} ? $options{-custom_plot_multiplier_5} : 10,
		'_custom_plot_data_5' => \%{$options{-custom_plot_data_5}},
		'_custom_plot_label_5' => exists $options{-custom_plot_label_5} ? $options{-custom_plot_label_5} : "Custom Plot 5",

		'_draw_custom_plot_6' => exists $options{-draw_custom_plot_6} ? $options{-draw_custom_plot_6} : 0,
		'_custom_plot_multiplier_6' => exists $options{-custom_plot_multiplier_6} ? $options{-custom_plot_multiplier_6} : 10,
		'_custom_plot_data_6' => \%{$options{-custom_plot_data_6}},
		'_custom_plot_label_6' => exists $options{-custom_plot_label_6} ? $options{-custom_plot_label_6} : "Custom Plot 6",

		'_draw_custom_plot_7' => exists $options{-draw_custom_plot_7} ? $options{-draw_custom_plot_7} : 0,
		'_custom_plot_multiplier_7' => exists $options{-custom_plot_multiplier_7} ? $options{-custom_plot_multiplier_7} : 10,
		'_custom_plot_data_7' => \%{$options{-custom_plot_data_7}},
		'_custom_plot_label_7' => exists $options{-custom_plot_label_7} ? $options{-custom_plot_label_7} : "Custom Plot 7",

		'_draw_custom_plot_8' => exists $options{-draw_custom_plot_8} ? $options{-draw_custom_plot_8} : 0,
		'_custom_plot_multiplier_8' => exists $options{-custom_plot_multiplier_8} ? $options{-custom_plot_multiplier_8} : 10,
		'_custom_plot_data_8' => \%{$options{-custom_plot_data_8}},
		'_custom_plot_label_8' => exists $options{-custom_plot_label_8} ? $options{-custom_plot_label_8} : "Custom Plot 8",

		'_topologies' => \%{$options{-topologies}},
		'_order' => \@{$options{-order}},
		'_bar_height' => exists $options{-bar_height} ? $options{-bar_height} : 15,
		'_seperator' => exists $options{-seperator} ? $options{-seperator} : 10,
		'_increment' => exists $options{-increment} ? $options{-increment} : 0,
		'_top' => exists $options{-top} ? $options{-top} : 40,
		'_multiplier' => exists $options{-multiplier} ? $options{-multiplier} : 1.6,
		'_text_buffer' => exists $options{-text_buffer} ? $options{-text_buffer} : 20,
		'_flank' => exists $options{-flank} ? $options{-flank} : 120,
		'_ttf_font' => exists $options{-ttf_font} ? $options{-ttf_font} : 0,
		'_ttf_font_size' => exists $options{-ttf_font_size} ? $options{-ttf_font_size} : 8,
		'_text_vertical_offset' => exists $options{-text_vertical_offset} ? $options{-text_vertical_offset} : 0,
		'_text_horizontal_offset' => exists $options{-text_horizontal_offset} ? $options{-text_horizontal_offset} : 0,
		'_signal_overides' => exists $options{-signal_overides} ? $options{-signal_overides} : 1,
		'_exclude_from_consensus' => exists $options{-exclude_from_consensus} ? $options{-exclude_from_consensus} : '',
		'_dont_colour' => exists $options{-dont_colour} ? $options{-dont_colour} : '',
		'_inside_rgb' => exists $options{-inside_rgb} ? $options{-inside_rgb} : [255,255,255],
		'_outside_rgb' => exists $options{-outside_rgb} ? $options{-outside_rgb} : [255,180,0],
		'_signal_rgb' => exists $options{-signal_rgb} ? $options{-signal_rgb} : [243,48,130]

	};

  	bless ($self,$class);
	return $self;

}

sub add_topology {

	my $self = shift;
	my %options = @_;

	if ($options{'-method'}){


		if ($options{'-n_term'}){
			$self->{'_topologies'}{$options{'-method'}}{'n_term'} = $options{'-n_term'};
		}else{
			$self->{'_topologies'}{$options{'-method'}}{'n_term'} = 'in';
		}

		if ($options{'-signal'}){
			$self->{'_topologies'}{$options{'-method'}}{'signal'} = $options{'-signal'};
		}

		@{$self->{'_topologies'}{$options{'-method'}}{'topology_array'}} = @{$options{'-topology'}} if exists ${$options{'-topology'}}[0];
		@{$self->{'_topologies'}{$options{'-method'}}{'reentrant_array'}} = @{$options{'-reentrant'}} if exists ${$options{'-reentrant'}}[0];

		$self->{'_topologies'}{$options{'-method'}}{'topology_string'} = $options{'-topology_string'} if $options{'-topology_string'};

		$self->{'_topologies'}{$options{'-method'}}{'reentrant_string'} = $options{'-reentrant_string'} if $options{'-reentrant_string'};
	}

	return $self;
}

sub png {

	my $self = shift;

	## Check options are numeric
	foreach (sort {$a cmp $b} keys %{$self->{'_topologies'}}){

		if ($self->{'_topologies'}{$_}{'topology_string'}){
			$self->{'_topologies'}{$_}{'topology_string'} =~ s/\D\.//g;
			$self->{'_topologies'}{$_}{'topology_string'} =~ s/;/,/g;
			@{$self->{'_topologies'}{$_}{'topology_array'}} = split(/,/,$self->{'_topologies'}{$_}{'topology_string'});

			foreach my $a (@{$self->{'_topologies'}{$_}{'topology_array'}}){
				die "\nHelix positions for $_ must be numeric.\n\n" if $a =~ /\D+/;
			}

			@{$self->{'_topologies'}{$_}{'topology_array'}} = sort {$a <=> $b} @{$self->{'_topologies'}{$_}{'topology_array'}};

		}elsif(scalar $self->{'_topologies'}{$_}{'topology'}){
			@{$self->{'_topologies'}{$_}{'topology_array'}} = @{$self->{'_topologies'}{$_}{'topology'}};
		}

		## Re-entrant stuff

		if ($self->{'_topologies'}{$_}{'reentrant_string'}){
			$self->{'_topologies'}{$_}{'reentrant_string'} =~ s/\D\.//g;
			$self->{'_topologies'}{$_}{'reentrant_string'} =~ s/;/,/g;
			@{$self->{'_topologies'}{$_}{'reentrant_array'}} = split(/,/,$self->{'_topologies'}{$_}{'reentrant_string'});

			foreach my $a (@{$self->{'_topologies'}{$_}{'reentrant_array'}}){
				die "\nRe-entrant helix positions for $_ must be numeric.\n\n" if $a =~ /\D+/;
			}

			@{$self->{'_topologies'}{$_}{'reentrant_array'}} = sort {$a <=> $b} @{$self->{'_topologies'}{$_}{'reentrant_array'}};

		}elsif(scalar $self->{'_topologies'}{$_}{'reentrant'}){
			@{$self->{'_topologies'}{$_}{'reentrant_array'}} = @{$self->{'_topologies'}{$_}{'reentrant'}};
		}

		###################

		if ($self->{'_topologies'}{$_}{'signal'}){
			die "\nSignal peptide position for $_ must be numeric.\n\n" if $self->{'_topologies'}{$_}{'signal'} =~ /\D+/;
		}

		## n-terminal defaults to outside in it's not in,inside,out,outside
		$self->{'_topologies'}{$_}{'n_term'} = 'out' if (($self->{'_topologies'}{$_}{'n_term'} ne 'in')&&($self->{'_topologies'}{$_}{'n_term'} ne 'inside')&&($self->{'_topologies'}{$_}{'n_term'} ne 'out')||$self->{'_topologies'}{$_}{'n_term'} eq 'outside');
		$self->{'_topologies'}{$_}{'n_term'} = 'in' if $self->{'_topologies'}{$_}{'n_term'} eq 'inside';
		$self->{'_topologies'}{$_}{'n_term'} = 'out' if $self->{'_topologies'}{$_}{'signal'};

	}

	my @numeric = ('_seq_length','_kd_window','_kd_multiplier','_bar_height','_seperator','_increment','_top','_multiplier','_text_buffer','_flank','_ttf_font_size','_text_vertical_offset','_text_horizontal_offset');

	foreach (@numeric){
		unless  ($self->{$_} =~ /[+-]?\d+\.\d+/){
			die "\nParameter $_ must be numeric.\n\n" if $self->{$_} =~ /[+-]?\D+/;
		}
	}

	if (exists $self->{'_ttf_font'}){

		$self->{'_ttf_font'} = 0 unless -e $self->{'_ttf_font'};

	}

	## Set a few things depending on what was provided
	if ($self->{'_method'} && (scalar $self->{'_topology_array'}||$self->{'_topology_string'})){

		$self->{'_topologies'}{$self->{'_method'}}{'n_term'} = $self->{'_n_term'};
		$self->{'_topologies'}{$self->{'_method'}}{'signal'} = $self->{'_signal'};
		$self->{'_topologies'}{$self->{'_method'}}{'topology_array'} = $self->{'_topology_array'} if scalar $self->{'_topology_array'};
		$self->{'_topologies'}{$self->{'_method'}}{'topology_string'} = $self->{'_topology_string'} if $self->{'_topology_string'};
	}

	unless ($self->{'_sequence'}){

		$self->{'_draw_kd'} = 0;

		if (!$self->{'_seq_length'}){
			$self->{'_seq_length'} = 500;
		}

	}else{
		$self->{'_sequence'} =~ s/\s+//g;
		$self->{'_seq_length'} = length $self->{'_sequence'};
	}

	$self->{'_draw_consensus'} = 0 if (scalar keys %{$self->{'_topologies'}} < 2);

	unless (scalar @{$self->{'_order'}}){

		@{$self->{'_order'}} = sort {$a cmp $b} keys %{$self->{'_topologies'}};

	}

	## Create image
	my $width = ($self->{'_seq_length'} + (1.5 * $self->{'_flank'})) * $self->{'_multiplier'};
	$width = 950 if $width < 950;

	my $height = (scalar @{$self->{'_order'}} * ($self->{'_bar_height'} + $self->{'_seperator'})) + ($self->{'_flank'}) - $self->{'_seperator'};
	$height = $height + 60 if $self->{'_draw_kd'};
	$height = $height + 60 if $self->{'_draw_custom_plot_1'};
	$height = $height + 60 if $self->{'_draw_custom_plot_2'};
	$height = $height + 60 if $self->{'_draw_custom_plot_3'};
	$height = $height + 60 if $self->{'_draw_custom_plot_4'};
	$height = $height + 60 if $self->{'_draw_custom_plot_5'};
	$height = $height + 60 if $self->{'_draw_custom_plot_6'};
	$height = $height + 60 if $self->{'_draw_custom_plot_7'};
	$height = $height + 60 if $self->{'_draw_custom_plot_8'};

	$height = $height + $self->{'_bar_height'} + $self->{'_seperator'}  if $self->{'_draw_consensus'};
	$self->{'im'} = new GD::Image($width,$height);

	## Check rgb colours are ok
	foreach (@{$self->{'_inside_rgb'}}){
		if ($_ =~ /\D+/ || $_ < 0 || $_ > 255){
			$self->{'_inside_rgb'} = [255,255,255];
		}
	}
	foreach (@{$self->{'_outside_rgb'}}){
		if ($_ =~ /\D+/ || $_ < 0 || $_ > 255){
			$self->{'_outside_rgb'} = [255,180,0];
		}
	}
	foreach (@{$self->{'_signal_rgb'}}){
		if ($_ =~ /\D+/ || $_ < 0 || $_ > 255){
			$self->{'_signal_rgb'} = [243,48,130];
		}
	}

	## Set up colours & some variables
	$self->{'inside'} = $self->{'im'}->colorAllocate(${$self->{'_inside_rgb'}}[0],${$self->{'_inside_rgb'}}[1],${$self->{'_inside_rgb'}}[2]);
	$self->{'outside'} = $self->{'im'}->colorAllocate(${$self->{'_outside_rgb'}}[0],${$self->{'_outside_rgb'}}[1],${$self->{'_outside_rgb'}}[2]);
	$self->{'signal'} = $self->{'im'}->colorAllocate(${$self->{'_signal_rgb'}}[0],${$self->{'_signal_rgb'}}[1],${$self->{'_signal_rgb'}}[2]);

	$self->{'black'} = $self->{'im'}->colorAllocate(0,0,0);
	$self->{'white'} = $self->{'im'}->colorAllocate(255,255,255);
	$self->{'dark_grey'}  = $self->{'im'}->colorAllocate(40,40,40);
	$self->{'dark_grey1'} = $self->{'im'}->colorAllocate(50,50,50);
	$self->{'dark_grey2'} = $self->{'im'}->colorAllocate(60,60,60);
	$self->{'dark_grey3'} = $self->{'im'}->colorAllocate(70,70,70);
	$self->{'green'} = $self->{'im'}->colorAllocate(37,188,25);

	$self->{'grey'} = $self->{'im'}->colorAllocate(190,190,190);
	$self->{'red'} = $self->{'im'}->colorAllocate(255,0,0);
	$self->{'x_start'} = 0;
	$self->{'y_start'} = 0;
	$self->{'x_stop'} = 0;
	$self->{'y_stop'} = 0;

	$self->{'im'}->fill(0,0,$self->{'white'});
	$self->write_text($self->{'_title'},6,4) if $self->{'_title'};

	## Main subroutines
	$self->draw_topologies;
	$self->draw_kd_plot if $self->{'_draw_kd'};
	$self->draw_custom_plot_1 if $self->{'_draw_custom_plot_1'};
	$self->draw_custom_plot_2 if $self->{'_draw_custom_plot_2'};
	$self->draw_custom_plot_3 if $self->{'_draw_custom_plot_3'};
	$self->draw_custom_plot_4 if $self->{'_draw_custom_plot_4'};
	$self->draw_custom_plot_5 if $self->{'_draw_custom_plot_5'};
	$self->draw_custom_plot_6 if $self->{'_draw_custom_plot_6'};
	$self->draw_custom_plot_7 if $self->{'_draw_custom_plot_7'};
	$self->draw_custom_plot_8 if $self->{'_draw_custom_plot_8'};
	$self->draw_consensus if $self->{'_draw_consensus'};
	$self->draw_grid if $self->{'_draw_grid'};
	$self->draw_key if $self->{'_draw_key'};

	return $self->{'im'}->GD::Image::png;

}

sub write_text {

	my $self = shift;
	my ($text,$x,$y) = @_;

	$x += $self->{'_text_horizontal_offset'};
	$y += $self->{'_text_vertical_offset'};

	if ($self->{'_ttf_font'}){
		$self->{'im'}->stringFT($self->{'black'},$self->{'_ttf_font'},$self->{'_ttf_font_size'},0,$x,$y + 12,$text,{linespacing=>0.6,charmap  => 'Unicode',});
	}else{
		$self->{'im'}->string(gdSmallFont,$x,$y,$text,$self->{'black'});
	}
}

sub draw_topologies {

	my $self = shift;

	#print "Drawing topology\n";

	my @dont_colour;
	if ($self->{'_dont_colour'}){
		@dont_colour = split(/,/,$self->{'_dont_colour'});
	}

	foreach (@{$self->{'_order'}}){



		my $dont_colour = 0;
		foreach my $dc (@dont_colour){
			$dont_colour++ if $_ eq $dc;
		}

		my $helix_count = 0;
		if (scalar @{$self->{'_topologies'}{$_}{'topology_array'}}){
			$helix_count = (scalar @{$self->{'_topologies'}{$_}{'topology_array'}}) / 2;
		}

		my $loop_count = $helix_count + 1;
		$loop_count = 0 if $helix_count == 0;
		$self->{'x_start'} = $self->{'_flank'};
		$self->{'y_start'} = $self->{'_top'} + $self->{'_increment'};
		$self->{'x_stop'} = $self->{'_flank'} + $self->{'_seq_length'};
		$self->{'y_stop'} = $self->{'_top'} + $self->{'_bar_height'} + $self->{'_increment'};

		$self->write_text($_,$self->{'_text_buffer'},$self->{'y_start'});

		$self->{'im'}->filledRectangle($self->{'x_start'} * $self->{'_multiplier'},$self->{'y_start'},$self->{'x_stop'} * $self->{'_multiplier'},$self->{'y_stop'},$self->{'black'});
		$self->{'im'}->filledRectangle(($self->{'x_start'} * $self->{'_multiplier'}) + 1,$self->{'y_start'} + 1,($self->{'x_stop'} * $self->{'_multiplier'}) - 1,$self->{'y_stop'} - 1,$self->{'inside'});

		my @tmp = @{$self->{'_topologies'}{$_}{'topology_array'}};
		my ($x1,$x2);
		my $count = 1;

	#print "topology: @tmp\n";

		## Draw Helices
		for (my $i = 0; $i < scalar @tmp;$i = $i+2){

			$tmp[$i] = 1 if $tmp[$i] < 1;
			$tmp[$i+1] = 1 if $tmp[$i] < 1;

			$x1 = ($tmp[$i] + $self->{'_flank'}) * $self->{'_multiplier'};
			$x2 = ($tmp[$i+1] + $self->{'_flank'}) * $self->{'_multiplier'};
			my $y1 = $self->{'y_start'};
			my $y2 = $self->{'y_stop'};

			$self->{'im'}->filledRectangle($x1,$y1,$x2,$y2,$self->{'black'});
			$x1++;
			$x2--;
			$y1++;
			$y2--;
			$self->{'im'}->filledRectangle($x1,$y1,$x2,$y2,$self->{'dark_grey3'});

		}

		## Draw Re-entrant Helices
		if (scalar \@{$self->{'_topologies'}{$_}{'reentrant_array'}}){

			my @re = @{$self->{'_topologies'}{$_}{'reentrant_array'}};

			for (my $i = 0; $i < scalar @re;$i = $i+2){

				$re[$i] = 1 if $re[$i] < 1;
				$re[$i+1] = 1 if $re[$i] < 1;

				$x1 = ($re[$i] + $self->{'_flank'}) * $self->{'_multiplier'};
				$x2 = ($re[$i+1] + $self->{'_flank'}) * $self->{'_multiplier'};
				my $y1 = $self->{'y_start'};
				my $y2 = $self->{'y_stop'};

				$self->{'im'}->filledRectangle($x1,$y1,$x2,$y2,$self->{'black'});
				$x1++;
				$x2--;
				$y1++;
				$y2--;
				$self->{'im'}->filledRectangle($x1,$y1,$x2,$y2,$self->{'green'});
			}
		}

		if ($self->{'_topologies'}{$_}{'n_term'} eq 'in'){

			for (my $i = 0; $i < scalar @tmp;$i = $i+2){

				next unless $tmp[$i+1];

				$x1 = $tmp[$i+1] + $self->{'_flank'};

				if (defined $tmp[$i+2]){
					$x2 = $tmp[$i+2] + $self->{'_flank'};
				}else{
					$x2 = $self->{'_seq_length'} + $self->{'_flank'};
				}

				unless($dont_colour){
					$self->{'im'}->filledRectangle((($x1) * $self->{'_multiplier'}) + 1,$self->{'y_start'} + 1,(($x2) * $self->{'_multiplier'}) - 1,$self->{'y_stop'} - 1,$self->{'outside'}) unless ($i / 2) % 2;
				}else{
					$self->{'im'}->filledRectangle((($x1) * $self->{'_multiplier'}) + 1,$self->{'y_start'} + 1,(($x2) * $self->{'_multiplier'}) - 1,$self->{'y_stop'} - 1,$self->{'white'}) unless ($i / 2) % 2;
				}

				$count++;
			}

		}else{

			for (my $i = 0; $i < scalar @tmp;$i = $i+2){

				$x1 = $tmp[$i-1] + $self->{'_flank'};
				if ($tmp[$i-1] > $tmp[$i]){
						$x1 = 0 + $self->{'_flank'};
				}
				$x2 = $tmp[$i] + $self->{'_flank'};

				unless($dont_colour){
					$self->{'im'}->filledRectangle((($x1) * $self->{'_multiplier'}) + 1,$self->{'y_start'} + 1,(($x2) * $self->{'_multiplier'}) - 1,$self->{'y_stop'} - 1,$self->{'outside'}) unless ($i / 2) % 2;
				}else{
					$self->{'im'}->filledRectangle((($x1) * $self->{'_multiplier'}) + 1,$self->{'y_start'} + 1,(($x2) * $self->{'_multiplier'}) - 1,$self->{'y_stop'} - 1,$self->{'white'}) unless ($i / 2) % 2;
				}
				$count++;
			}

			if (defined $tmp[-1]){
				$x1 = $tmp[-1] + $self->{'_flank'};
			}else{
				$x1 = 0 + $self->{'_flank'};
			}

			unless($dont_colour){
				$self->{'im'}->filledRectangle(($x1 * $self->{'_multiplier'}) + 1,$self->{'y_start'} + 1,(($self->{'_seq_length'} + $self->{'_flank'}) * $self->{'_multiplier'}) - 1,$self->{'y_stop'} - 1,$self->{'outside'}) unless (scalar @tmp / 2) % 2;
			}else{
				$self->{'im'}->filledRectangle(($x1 * $self->{'_multiplier'}) + 1,$self->{'y_start'} + 1,(($self->{'_seq_length'} + $self->{'_flank'}) * $self->{'_multiplier'}) - 1,$self->{'y_stop'} - 1,$self->{'white'}) unless (scalar @tmp / 2) % 2;
			}
		}

		## Draw signal peptide
		if ($self->{'_topologies'}{$_}{'signal'}){
				$self->{'im'}->filledRectangle($self->{'_flank'} * $self->{'_multiplier'},$self->{'y_start'},(($self->{'_flank'} + $self->{'_topologies'}{$_}{'signal'}) * $self->{'_multiplier'}),$self->{'y_stop'},$self->{'black'});
				$self->{'im'}->filledRectangle(($self->{'_flank'} * $self->{'_multiplier'}) + 1,$self->{'y_start'} + 1,(($self->{'_flank'} + $self->{'_topologies'}{$_}{'signal'}) * $self->{'_multiplier'}) - 1,$self->{'y_stop'} - 1,$self->{'signal'});
		}

		## Draw Re-entrant Helices
		if (scalar \@{$self->{'_topologies'}{$_}{'reentrant_array'}}){

			my @re = @{$self->{'_topologies'}{$_}{'reentrant_array'}};

			for (my $i = 0; $i < scalar @re;$i = $i+2){

				$re[$i] = 1 if $re[$i] < 1;
				$re[$i+1] = 1 if $re[$i] < 1;

				$x1 = ($re[$i] + $self->{'_flank'}) * $self->{'_multiplier'};
				$x2 = ($re[$i+1] + $self->{'_flank'}) * $self->{'_multiplier'};
				my $y1 = $self->{'y_start'};
				my $y2 = $self->{'y_stop'};

				$self->{'im'}->filledRectangle($x1,$y1,$x2,$y2,$self->{'black'});
				$x1++;
				$x2--;
				$y1++;
				$y2--;
				$self->{'im'}->filledRectangle($x1,$y1,$x2,$y2,$self->{'green'});
			}
		}


		$self->{'_increment'} += $self->{'_bar_height'} + $self->{'_seperator'};
	}

}

sub draw_kd_plot {

	my $self = shift;

	$self->generate_kd_data;

	$self->{'_kd_multiplier'} = $self->{'_kd_multiplier'} * -1;

	$self->write_text($self->{'_kd_label'},$self->{'_text_buffer'},$self->{'_increment'} + $self->{'_top'} + 22);

	#$self->{'im'}->setThickness(1);
	$self->{'im'}->dashedLine($self->{'_flank'} * $self->{'_multiplier'},$self->{'_increment'} + $self->{'_top'} + 30,($self->{'_flank'} + $self->{'_seq_length'}) * $self->{'_multiplier'},$self->{'_increment'} + $self->{'_top'} + 30,$self->{'dark_grey3'});

	my @tmp = sort {$a<=>$b} keys %{$self->{'kd_data'}};

	foreach (sort {$a<=>$b} keys %{$self->{'kd_data'}}){

		my $x1 = ($self->{'_flank'} + $_) * $self->{'_multiplier'};
		my $y1 = ($self->{'_kd_multiplier'} * $self->{'kd_data'}{$_}) + $self->{'_increment'} + $self->{'_top'} + 30;
		my ($x2,$y2);

		if ($_ == $tmp[-1]){
			$x2 = $x1;
			$y2 = $y1;
		}else{
			$x2 = ($self->{'_flank'} + ($_ + 1)) * $self->{'_multiplier'};
			$y2 = ($self->{'_kd_multiplier'} * $self->{'kd_data'}{$_ + 1}) + $self->{'_increment'} + $self->{'_top'} + 30;
		}

		$self->{'im'}->line($x1,$y1,$x2,$y2,$self->{'dark_grey1'})
	}
}

sub draw_custom_plot_1 {

	my $self = shift;

	$self->{'_increment'} += 60 if $self->{'_draw_custom_plot_1'};

	$self->{'_custom_plot_multiplier_1'} = $self->{'_custom_plot_multiplier_1'} * -1;

	$self->write_text($self->{'_custom_plot_label_1'},$self->{'_text_buffer'},$self->{'_increment'} + $self->{'_top'} + 22);

	#$self->{'im'}->setThickness(1);
	$self->{'im'}->dashedLine($self->{'_flank'} * $self->{'_multiplier'},$self->{'_increment'} + $self->{'_top'} + 30,($self->{'_flank'} + $self->{'_seq_length'}) * $self->{'_multiplier'},$self->{'_increment'} + $self->{'_top'} + 30,$self->{'dark_grey3'});



	my @tmp = sort {$a<=>$b} keys %{$self->{'_custom_plot_data_1'}};

	foreach (sort {$a<=>$b} keys %{$self->{'_custom_plot_data_1'}}){

		my $x1 = ($self->{'_flank'} + $_) * $self->{'_multiplier'};
		my $y1 = ($self->{'_custom_plot_multiplier_1'} * $self->{'_custom_plot_data_1'}{$_}) + $self->{'_increment'} + $self->{'_top'} + 30;
		my ($x2,$y2);

		if ($_ == $tmp[-1]){
			$x2 = $x1;
			$y2 = $y1;
		}else{
			$x2 = ($self->{'_flank'} + ($_ + 1)) * $self->{'_multiplier'};
			$y2 = ($self->{'_custom_plot_multiplier_1'} * $self->{'_custom_plot_data_1'}{$_ + 1}) + $self->{'_increment'} + $self->{'_top'} + 30;
		}

		$self->{'im'}->line($x1,$y1,$x2,$y2,$self->{'dark_grey1'})
	}

	## Draw line at 0.5
	#my $point_five = ($self->{'_custom_plot_multiplier_1'} * 0.5) + $self->{'_increment'} + $self->{'_top'} + 30;
	#$self->{'im'}->dashedLine($self->{'_flank'} * $self->{'_multiplier'},$point_five,($self->{'_flank'} + $self->{'_seq_length'}) * $self->{'_multiplier'},$point_five,$self->{'red'});

}

sub draw_custom_plot_2 {

	my $self = shift;

	$self->{'_increment'} += 60 if $self->{'_draw_custom_plot_2'};

	$self->{'_custom_plot_multiplier_2'} = $self->{'_custom_plot_multiplier_2'} * -1;

	$self->write_text($self->{'_custom_plot_label_2'},$self->{'_text_buffer'},$self->{'_increment'} + $self->{'_top'} + 22);

	#$self->{'im'}->setThickness(1);
	$self->{'im'}->dashedLine($self->{'_flank'} * $self->{'_multiplier'},$self->{'_increment'} + $self->{'_top'} + 30,($self->{'_flank'} + $self->{'_seq_length'}) * $self->{'_multiplier'},$self->{'_increment'} + $self->{'_top'} + 30,$self->{'dark_grey3'});



	my @tmp = sort {$a<=>$b} keys %{$self->{'_custom_plot_data_2'}};

	foreach (sort {$a<=>$b} keys %{$self->{'_custom_plot_data_2'}}){

		my $x1 = ($self->{'_flank'} + $_) * $self->{'_multiplier'};
		my $y1 = ($self->{'_custom_plot_multiplier_2'} * $self->{'_custom_plot_data_2'}{$_}) + $self->{'_increment'} + $self->{'_top'} + 30;
		my ($x2,$y2);

		if ($_ == $tmp[-1]){
			$x2 = $x1;
			$y2 = $y1;
		}else{
			$x2 = ($self->{'_flank'} + ($_ + 1)) * $self->{'_multiplier'};
			$y2 = ($self->{'_custom_plot_multiplier_2'} * $self->{'_custom_plot_data_2'}{$_ + 1}) + $self->{'_increment'} + $self->{'_top'} + 30;
		}

		$self->{'im'}->line($x1,$y1,$x2,$y2,$self->{'dark_grey1'})
	}

	## Draw line at 0.211197
	my $point_five = ($self->{'_custom_plot_multiplier_2'} * 0.211197) + $self->{'_increment'} + $self->{'_top'} + 30;
	$self->{'im'}->dashedLine($self->{'_flank'} * $self->{'_multiplier'},$point_five,($self->{'_flank'} + $self->{'_seq_length'}) * $self->{'_multiplier'},$point_five,$self->{'red'});

}

sub draw_custom_plot_3 {

	my $self = shift;

	$self->{'_increment'} += 60 if $self->{'_draw_custom_plot_3'};

	$self->{'_custom_plot_multiplier_3'} = $self->{'_custom_plot_multiplier_3'} * -1;

	$self->write_text($self->{'_custom_plot_label_3'},$self->{'_text_buffer'},$self->{'_increment'} + $self->{'_top'} + 22);

	#$self->{'im'}->setThickness(1);
	$self->{'im'}->dashedLine($self->{'_flank'} * $self->{'_multiplier'},$self->{'_increment'} + $self->{'_top'} + 30,($self->{'_flank'} + $self->{'_seq_length'}) * $self->{'_multiplier'},$self->{'_increment'} + $self->{'_top'} + 30,$self->{'dark_grey3'});



	my @tmp = sort {$a<=>$b} keys %{$self->{'_custom_plot_data_3'}};

	foreach (sort {$a<=>$b} keys %{$self->{'_custom_plot_data_3'}}){

		my $x1 = ($self->{'_flank'} + $_) * $self->{'_multiplier'};
		my $y1 = ($self->{'_custom_plot_multiplier_3'} * $self->{'_custom_plot_data_3'}{$_}) + $self->{'_increment'} + $self->{'_top'} + 30;
		my ($x2,$y2);

		if ($_ == $tmp[-1]){
			$x2 = $x1;
			$y2 = $y1;
		}else{
			$x2 = ($self->{'_flank'} + ($_ + 1)) * $self->{'_multiplier'};
			$y2 = ($self->{'_custom_plot_multiplier_3'} * $self->{'_custom_plot_data_3'}{$_ + 1}) + $self->{'_increment'} + $self->{'_top'} + 30;
		}

		$self->{'im'}->line($x1,$y1,$x2,$y2,$self->{'dark_grey1'})
	}

	## Draw line at 0.5
	#my $point_five = ($self->{'_custom_plot_multiplier_3'} * 0.5) + $self->{'_increment'} + $self->{'_top'} + 30;
	#$self->{'im'}->dashedLine($self->{'_flank'} * $self->{'_multiplier'},$point_five,($self->{'_flank'} + $self->{'_seq_length'}) * $self->{'_multiplier'},$point_five,$self->{'red'});

}

sub draw_custom_plot_4 {

	my $self = shift;

	$self->{'_increment'} += 60 if $self->{'_draw_custom_plot_4'};

	$self->{'_custom_plot_multiplier_4'} = $self->{'_custom_plot_multiplier_4'} * -1;

	$self->write_text($self->{'_custom_plot_label_4'},$self->{'_text_buffer'},$self->{'_increment'} + $self->{'_top'} + 22);

	#$self->{'im'}->setThickness(1);
	$self->{'im'}->dashedLine($self->{'_flank'} * $self->{'_multiplier'},$self->{'_increment'} + $self->{'_top'} + 30,($self->{'_flank'} + $self->{'_seq_length'}) * $self->{'_multiplier'},$self->{'_increment'} + $self->{'_top'} + 30,$self->{'dark_grey3'});



	my @tmp = sort {$a<=>$b} keys %{$self->{'_custom_plot_data_4'}};

	foreach (sort {$a<=>$b} keys %{$self->{'_custom_plot_data_4'}}){

		my $x1 = ($self->{'_flank'} + $_) * $self->{'_multiplier'};
		my $y1 = ($self->{'_custom_plot_multiplier_4'} * $self->{'_custom_plot_data_4'}{$_}) + $self->{'_increment'} + $self->{'_top'} + 30;
		my ($x2,$y2);

		if ($_ == $tmp[-1]){
			$x2 = $x1;
			$y2 = $y1;
		}else{
			$x2 = ($self->{'_flank'} + ($_ + 1)) * $self->{'_multiplier'};
			$y2 = ($self->{'_custom_plot_multiplier_4'} * $self->{'_custom_plot_data_4'}{$_ + 1}) + $self->{'_increment'} + $self->{'_top'} + 30;
		}

		$self->{'im'}->line($x1,$y1,$x2,$y2,$self->{'dark_grey1'})
	}

	## Draw line at 0.5
	my $point_five = ($self->{'_custom_plot_multiplier_4'} * 0.5) + $self->{'_increment'} + $self->{'_top'} + 30;
	$self->{'im'}->dashedLine($self->{'_flank'} * $self->{'_multiplier'},$point_five,($self->{'_flank'} + $self->{'_seq_length'}) * $self->{'_multiplier'},$point_five,$self->{'red'});

}

sub draw_custom_plot_5 {

	my $self = shift;

	$self->{'_increment'} += 60 if $self->{'_draw_custom_plot_5'};

	$self->{'_custom_plot_multiplier_5'} = $self->{'_custom_plot_multiplier_5'} * -1;

	$self->write_text($self->{'_custom_plot_label_5'},$self->{'_text_buffer'},$self->{'_increment'} + $self->{'_top'} + 22);

	#$self->{'im'}->setThickness(1);
	$self->{'im'}->dashedLine($self->{'_flank'} * $self->{'_multiplier'},$self->{'_increment'} + $self->{'_top'} + 50,($self->{'_flank'} + $self->{'_seq_length'}) * $self->{'_multiplier'},$self->{'_increment'} + $self->{'_top'} + 50,$self->{'dark_grey3'});



	my @tmp = sort {$a<=>$b} keys %{$self->{'_custom_plot_data_5'}};

	foreach (sort {$a<=>$b} keys %{$self->{'_custom_plot_data_5'}}){

		my $x1 = ($self->{'_flank'} + $_) * $self->{'_multiplier'};
		my $y1 = ($self->{'_custom_plot_multiplier_5'} * $self->{'_custom_plot_data_5'}{$_}) + $self->{'_increment'} + $self->{'_top'} + 50;
		my ($x2,$y2);

		if ($_ == $tmp[-1]){
			$x2 = $x1;
			$y2 = $y1;
		}else{
			$x2 = ($self->{'_flank'} + ($_ + 1)) * $self->{'_multiplier'};
			$y2 = ($self->{'_custom_plot_multiplier_5'} * $self->{'_custom_plot_data_5'}{$_ + 1}) + $self->{'_increment'} + $self->{'_top'} + 50;
		}

		$self->{'im'}->line($x1,$y1,$x2,$y2,$self->{'dark_grey1'})
	}

	## Draw line at 0.5
	#my $point_five = ($self->{'_custom_plot_multiplier_5'} * 0.5) + $self->{'_increment'} + $self->{'_top'} + 50;
	#$self->{'im'}->dashedLine($self->{'_flank'} * $self->{'_multiplier'},$point_five,($self->{'_flank'} + $self->{'_seq_length'}) * $self->{'_multiplier'},$point_five,$self->{'red'});

}

sub draw_custom_plot_6 {

	my $self = shift;

	$self->{'_increment'} += 60 if $self->{'_draw_custom_plot_6'};

	$self->{'_custom_plot_multiplier_6'} = $self->{'_custom_plot_multiplier_6'} * -1;

	$self->write_text($self->{'_custom_plot_label_6'},$self->{'_text_buffer'},$self->{'_increment'} + $self->{'_top'} + 22);

	#$self->{'im'}->setThickness(1);
	$self->{'im'}->dashedLine($self->{'_flank'} * $self->{'_multiplier'},$self->{'_increment'} + $self->{'_top'} + 30,($self->{'_flank'} + $self->{'_seq_length'}) * $self->{'_multiplier'},$self->{'_increment'} + $self->{'_top'} + 30,$self->{'dark_grey3'});



	my @tmp = sort {$a<=>$b} keys %{$self->{'_custom_plot_data_6'}};

	foreach (sort {$a<=>$b} keys %{$self->{'_custom_plot_data_6'}}){

		my $x1 = ($self->{'_flank'} + $_) * $self->{'_multiplier'};
		my $y1 = ($self->{'_custom_plot_multiplier_6'} * $self->{'_custom_plot_data_6'}{$_}) + $self->{'_increment'} + $self->{'_top'} + 30;
		my ($x2,$y2);

		if ($_ == $tmp[-1]){
			$x2 = $x1;
			$y2 = $y1;
		}else{
			$x2 = ($self->{'_flank'} + ($_ + 1)) * $self->{'_multiplier'};
			$y2 = ($self->{'_custom_plot_multiplier_6'} * $self->{'_custom_plot_data_6'}{$_ + 1}) + $self->{'_increment'} + $self->{'_top'} + 30;
		}

		$self->{'im'}->line($x1,$y1,$x2,$y2,$self->{'dark_grey1'})
	}

	## Draw line at 0.5
	my $point_five = ($self->{'_custom_plot_multiplier_6'} * 0.5) + $self->{'_increment'} + $self->{'_top'} + 30;
	$self->{'im'}->dashedLine($self->{'_flank'} * $self->{'_multiplier'},$point_five,($self->{'_flank'} + $self->{'_seq_length'}) * $self->{'_multiplier'},$point_five,$self->{'red'});

}

sub draw_custom_plot_7 {

	my $self = shift;

	$self->{'_increment'} += 60 if $self->{'_draw_custom_plot_7'};

	$self->{'_custom_plot_multiplier_7'} = $self->{'_custom_plot_multiplier_7'} * -1;

	$self->write_text($self->{'_custom_plot_label_7'},$self->{'_text_buffer'},$self->{'_increment'} + $self->{'_top'} + 22);

	#$self->{'im'}->setThickness(1);
	$self->{'im'}->dashedLine($self->{'_flank'} * $self->{'_multiplier'},$self->{'_increment'} + $self->{'_top'} + 50,($self->{'_flank'} + $self->{'_seq_length'}) * $self->{'_multiplier'},$self->{'_increment'} + $self->{'_top'} + 50,$self->{'dark_grey3'});



	my @tmp = sort {$a<=>$b} keys %{$self->{'_custom_plot_data_7'}};

	foreach (sort {$a<=>$b} keys %{$self->{'_custom_plot_data_7'}}){

		my $x1 = ($self->{'_flank'} + $_) * $self->{'_multiplier'};
		my $y1 = ($self->{'_custom_plot_multiplier_7'} * $self->{'_custom_plot_data_7'}{$_}) + $self->{'_increment'} + $self->{'_top'} + 50;
		my ($x2,$y2);

		if ($_ == $tmp[-1]){
			$x2 = $x1;
			$y2 = $y1;
		}else{
			$x2 = ($self->{'_flank'} + ($_ + 1)) * $self->{'_multiplier'};
			$y2 = ($self->{'_custom_plot_multiplier_7'} * $self->{'_custom_plot_data_7'}{$_ + 1}) + $self->{'_increment'} + $self->{'_top'} + 50;
		}

		$self->{'im'}->line($x1,$y1,$x2,$y2,$self->{'dark_grey1'})
	}

	## Draw line at 0.5
	#my $point_five = ($self->{'_custom_plot_multiplier_5'} * 0.5) + $self->{'_increment'} + $self->{'_top'} + 50;
	#$self->{'im'}->dashedLine($self->{'_flank'} * $self->{'_multiplier'},$point_five,($self->{'_flank'} + $self->{'_seq_length'}) * $self->{'_multiplier'},$point_five,$self->{'red'});

}

sub draw_custom_plot_8 {

	my $self = shift;

	$self->{'_increment'} += 60 if $self->{'_draw_custom_plot_8'};

	$self->{'_custom_plot_multiplier_8'} = $self->{'_custom_plot_multiplier_8'} * -1;

	$self->write_text($self->{'_custom_plot_label_8'},$self->{'_text_buffer'},$self->{'_increment'} + $self->{'_top'} + 22);

	#$self->{'im'}->setThickness(1);
	$self->{'im'}->dashedLine($self->{'_flank'} * $self->{'_multiplier'},$self->{'_increment'} + $self->{'_top'} + 30,($self->{'_flank'} + $self->{'_seq_length'}) * $self->{'_multiplier'},$self->{'_increment'} + $self->{'_top'} + 30,$self->{'dark_grey3'});



	my @tmp = sort {$a<=>$b} keys %{$self->{'_custom_plot_data_8'}};

	foreach (sort {$a<=>$b} keys %{$self->{'_custom_plot_data_8'}}){

		my $x1 = ($self->{'_flank'} + $_) * $self->{'_multiplier'};
		my $y1 = ($self->{'_custom_plot_multiplier_8'} * $self->{'_custom_plot_data_8'}{$_}) + $self->{'_increment'} + $self->{'_top'} + 30;
		my ($x2,$y2);

		if ($_ == $tmp[-1]){
			$x2 = $x1;
			$y2 = $y1;
		}else{
			$x2 = ($self->{'_flank'} + ($_ + 1)) * $self->{'_multiplier'};
			$y2 = ($self->{'_custom_plot_multiplier_8'} * $self->{'_custom_plot_data_8'}{$_ + 1}) + $self->{'_increment'} + $self->{'_top'} + 30;
		}

		$self->{'im'}->line($x1,$y1,$x2,$y2,$self->{'dark_grey1'})
	}

	## Draw line at 0.5
	my $point_five = ($self->{'_custom_plot_multiplier_8'} * 0.5) + $self->{'_increment'} + $self->{'_top'} + 30;
	$self->{'im'}->dashedLine($self->{'_flank'} * $self->{'_multiplier'},$point_five,($self->{'_flank'} + $self->{'_seq_length'}) * $self->{'_multiplier'},$point_five,$self->{'red'});

}

sub generate_kd_data{

	my $self = shift;

	# KD index
	my %kd = (
	'A' =>  1.800,
	'R' => -4.500,
	'N' => -3.500,
	'D' => -3.500,
	'C' =>  2.500,
	'Q' => -3.500,
	'E' => -3.500,
	'G' => -0.400,
	'H' => -3.200,
	'I' =>  4.500,
	'L' =>  3.800,
	'K' => -3.900,
	'M' =>  1.900,
	'F' =>  2.800,
	'P' => -1.600,
	'S' => -0.800,
	'T' => -0.700,
	'W' => -0.900,
	'Y' => -1.300,
	'V' =>  4.200
	);

	my %hash;
	$self->{'_kd_window'}++ unless $self->{'_kd_window'} % 2;
	$self->{'_sequence'} =~ s/\s+//g;
	my $start = ($self->{'_kd_window'} / 2) - 0.5;
	my $limit = $self->{'_seq_length'} - $start;

	for (my $i = $start; $i < $limit; $i++){

		my $ave = 0;

		for (my $j = ($i - $start); $j < ($i + $start + 1); $j++){

			my $s = uc substr($self->{'_sequence'},$j,1);
			$ave += $kd{$s} if exists $kd{$s};

		}

		$hash{$i+1} = $ave/$self->{'_kd_window'};
	}

	%{$self->{'kd_data'}} = %hash;

}

sub draw_consensus {

	my $self = shift;

	my $consensus_hash = ();
	for (1..$self->{'_seq_length'}){
		$$consensus_hash{$_}{'h'} = 0;
		$$consensus_hash{$_}{'i'} = 0;
		$$consensus_hash{$_}{'o'} = 0;
		$$consensus_hash{$_}{'s'} = 0;
	}

	my @exclude;
	if ($self->{'_exclude_from_consensus'}){
		@exclude = split(/,/,$self->{'_exclude_from_consensus'});
	}

	for my $p (1..$self->{'_seq_length'}){

		my $type = '';

		foreach my $t (keys %{$self->{'_topologies'}}){

			my @tmp = @{$self->{'_topologies'}{$t}{'topology_array'}};

			my $count = 0;

			for (my $i = 0; $i < scalar @tmp;$i = $i+2){

				$count++;


				if ($p >= $tmp[$i] && $p <= $tmp[$i+1]){

					$type = 'helix';

				}else{

					my ($previous,$next);

					if (defined $tmp[$i-1]){
						$previous = $tmp[$i-1];
						$previous = 0 if $previous > $tmp[$i];
					}else{
						$previous = 0;
					}

					if (defined $tmp[$i+2]){
						$next = $tmp[$i+2];
						$next = ($self->{'_seq_length'} + 1) if $next < $tmp[$i];
					}else{
						$next = $self->{'_seq_length'} + 1;
					}

					if ($count % 2){

						$type = 'in'  if $p < $tmp[$i] && $p > $previous && $self->{'_topologies'}{$t}{'n_term'} eq 'in';
						$type = 'out'  if $p > $tmp[$i+1] && $p < $next && $self->{'_topologies'}{$t}{'n_term'} eq 'in';
						$type = 'out'  if $p < $tmp[$i] && $p > $previous && $self->{'_topologies'}{$t}{'n_term'} eq 'out';
						$type = 'in'  if $p > $tmp[$i+1] && $p < $next && $self->{'_topologies'}{$t}{'n_term'} eq 'out';

					}else{

						$type = 'out'  if $p < $tmp[$i] && $p > $previous && $self->{'_topologies'}{$t}{'n_term'} eq 'in';
						$type = 'in'  if $p > $tmp[$i+1] && $p < $next && $self->{'_topologies'}{$t}{'n_term'} eq 'in';
						$type = 'in'  if $p < $tmp[$i] && $p > $previous && $self->{'_topologies'}{$t}{'n_term'} eq 'out';
						$type = 'out'  if $p > $tmp[$i+1] && $p < $next && $self->{'_topologies'}{$t}{'n_term'} eq 'out';
					}
				}
			}

			$type = $self->{'_topologies'}{$t}{'n_term'} if scalar @tmp == 0;

			my $skip = 0;

			if ($self->{'_topologies'}{$t}{'signal'}){
				if (($p <= $self->{'_topologies'}{$t}{'signal'}) && ($self->{'_signal_overides'})){
					$$consensus_hash{$p}{'s'}++;
					$skip++;
				}
			}
			next if $skip;

			foreach my $e (@exclude){
				$skip++ if $e eq $t;
			}
			next if $skip;

			$$consensus_hash{$p}{'h'}++ if $type eq 'helix';
			$$consensus_hash{$p}{'i'}++ if $type eq 'in';
			$$consensus_hash{$p}{'o'}++ if $type eq 'out';


		}
	}

	$self->{'_increment'} += 60 if $self->{'_draw_kd'};

	$self->{'x_start'} = $self->{'_flank'};
	$self->{'y_start'} = $self->{'_top'} + $self->{'_increment'};
	$self->{'x_stop'} = $self->{'_flank'} + $self->{'_seq_length'};
	$self->{'y_stop'} = $self->{'_top'} + $self->{'_bar_height'} + $self->{'_increment'};

	$self->write_text($self->{'_consensus_label'},$self->{'_text_buffer'},$self->{'y_start'});

	$self->{'im'}->filledRectangle($self->{'x_start'} * $self->{'_multiplier'},$self->{'y_start'},$self->{'x_stop'} * $self->{'_multiplier'},$self->{'y_stop'},$self->{'black'});
	$self->{'im'}->filledRectangle(($self->{'x_start'} * $self->{'_multiplier'}) + 1,$self->{'y_start'} + 1,($self->{'x_stop'} * $self->{'_multiplier'}) - 1,$self->{'y_stop'} - 1,$self->{'inside'});


	my $consensus_string;
	foreach (sort {$a<=>$b} keys  %{$consensus_hash}){

		# print consensus details
		# printf "$_ h: %d i: %d o: %d s: %d\n",$$consensus_hash{$_}{'h'},$$consensus_hash{$_}{'i'},$$consensus_hash{$_}{'o'},$$consensus_hash{$_}{'s'};

		if ($$consensus_hash{$_}{'h'} > $$consensus_hash{$_}{'i'} && $$consensus_hash{$_}{'h'} > $$consensus_hash{$_}{'o'} && $$consensus_hash{$_}{'h'} > $$consensus_hash{$_}{'s'}){
			$consensus_string .= 'h';
		}elsif($$consensus_hash{$_}{'i'} > $$consensus_hash{$_}{'o'} && $$consensus_hash{$_}{'i'} > $$consensus_hash{$_}{'h'} && $$consensus_hash{$_}{'i'} > $$consensus_hash{$_}{'s'}){
			$consensus_string .= 'i';
		}elsif($$consensus_hash{$_}{'o'} > $$consensus_hash{$_}{'h'} && $$consensus_hash{$_}{'o'} > $$consensus_hash{$_}{'i'} && $$consensus_hash{$_}{'o'} > $$consensus_hash{$_}{'s'}){
			$consensus_string .= 'o';
		}elsif($$consensus_hash{$_}{'s'} > $$consensus_hash{$_}{'h'} && $$consensus_hash{$_}{'s'} > $$consensus_hash{$_}{'i'} && $$consensus_hash{$_}{'s'} > $$consensus_hash{$_}{'o'}){
			$consensus_string .= 's';
		}else{
			$consensus_string .= '-';
		}
	}

	my @array = split (//,$consensus_string);
	my $count = 1;
	my @tmp;
	my $in_helix = 0;
	my $consensus_helix = 0;
	foreach (@array){

		if (($_ eq 'h')&&($in_helix == 0)){
			$in_helix = 1;
			push @tmp,$count;
			$consensus_helix++;
		}
		if (($_ ne 'h')&&($in_helix == 1)){
			$in_helix = 0;
			push @tmp,($count - 1);
			$consensus_helix++;
		}
		$count++;
	}

	## Draw loops
	$count = 1;
	foreach my $p (@array){

		my $x1 = ($count - 1 + $self->{'_flank'}) * $self->{'_multiplier'};
		my $x2 = ($count + 1 + $self->{'_flank'}) * $self->{'_multiplier'};
		my $y1 = $self->{'y_start'};
		my $y2 = $self->{'y_stop'};
		$x1++;
		$x2--;
		$y1++;
		$y2--;

		if ($count == $self->{'_seq_length'}){
			$x2 = $x2 - 2;
		}

		if ($p eq 'i'){
			$self->{'im'}->filledRectangle($x1,$y1,$x2,$y2,$self->{'inside'});
		}elsif ($p eq 'o'){
			$self->{'im'}->filledRectangle($x1,$y1,$x2,$y2,$self->{'outside'});
		}elsif ($p eq 's'){
			$self->{'im'}->filledRectangle($x1,$y1,$x2,$y2,$self->{'signal'});
		}elsif ($p eq '-'){
			$self->{'im'}->filledRectangle($x1,$y1,$x2,$y2,$self->{'grey'});
		}
		$count++;
	}

	## Draw Helices
	for (my $i = 0; $i < (scalar @tmp);$i = $i+2){

		my $x1 = ($tmp[$i] + $self->{'_flank'}) * $self->{'_multiplier'};
		my $x2 = ($tmp[$i+1] + $self->{'_flank'}) * $self->{'_multiplier'};
		my $y1 = $self->{'y_start'};
		my $y2 = $self->{'y_stop'};

		$self->{'im'}->filledRectangle($x1,$y1,$x2,$y2,$self->{'black'});
		$x1++;
		$x2--;
		$y1++;
		$y2--;
		$self->{'im'}->filledRectangle($x1,$y1,$x2,$y2,$self->{'dark_grey3'});

	}
}

sub draw_grid {

	my $self = shift;
	my ($n,$m) = split(/\./,($self->{'_seq_length'}/100));

	for (1..$n){

		my $x1 = ($self->{'_flank'} + ($_*100)) * $self->{'_multiplier'};
		my $y1 = $self->{'_top'} - 10;

		my $x2 = ($self->{'_flank'} + ($_*100)) * $self->{'_multiplier'};
		my $y2 = ((scalar (@{$self->{'_order'}}) * ($self->{'_bar_height'} + $self->{'_seperator'})) + $self->{'_top'});

		$y2 = $y2 + $self->{'_bar_height'} + $self->{'_seperator'}  if $self->{'_draw_consensus'};
		$y2 = $y2 + 60 if $self->{'_draw_kd'};
		$y2 = $y2 + 60 if $self->{'_draw_custom_plot'};

		$self->{'im'}->line($x1,$y1,$x2,$y2,$self->{'red'});

		my $string = $_ * 100;

		$x1 = $x1 - 8;
		$y1 = $y1 - 14;

		$self->write_text($string,$x1,$y1);
	}
}

sub draw_key {

	my $self = shift;
	my $key_box = 15;
	my $mod = 0;

	my $y1 = ((scalar (@{$self->{'_order'}}) * ($self->{'_bar_height'} + $self->{'_seperator'})) + $self->{'_top'});
	$y1 = $y1 + $self->{'_bar_height'} + $self->{'_seperator'}  if $self->{'_draw_consensus'};
	$y1 = $y1 + 60 if $self->{'_draw_kd'};

	$y1 = $y1 + 60 if $self->{'_draw_custom_plot_1'};
	$y1 = $y1 + 20 if $self->{'_draw_custom_plot_1'};
	$y1 = $y1 + 60 if $self->{'_draw_custom_plot_2'};
	$y1 = $y1 + 60 if $self->{'_draw_custom_plot_3'};
	$y1 = $y1 + 60 if $self->{'_draw_custom_plot_4'};
	$y1 = $y1 + 60 if $self->{'_draw_custom_plot_5'};
	$y1 = $y1 + 60 if $self->{'_draw_custom_plot_6'};
	$y1 = $y1 + 60 if $self->{'_draw_custom_plot_7'};
	$y1 = $y1 + 60 if $self->{'_draw_custom_plot_8'};

	$self->{'im'}->filledRectangle(($self->{'_flank'} * $self->{'_multiplier'}) + $mod,$y1 + 20,($self->{'_flank'} * $self->{'_multiplier'}) + $key_box + $mod,$y1 + 20 + $key_box,$self->{'black'});
	$self->{'im'}->filledRectangle(($self->{'_flank'} * $self->{'_multiplier'}) + 1 + $mod,$y1 + 21,($self->{'_flank'} * $self->{'_multiplier'}) + $key_box - 1 + $mod,$y1 + 20 + $key_box - 1,$self->{'signal'});
	$self->write_text($self->{'_signal_label'},($self->{'_flank'} * $self->{'_multiplier'}) + $key_box + 5 + $mod,$y1 + 6 + $key_box);

	$mod += 150;

	$self->{'im'}->filledRectangle(($self->{'_flank'} * $self->{'_multiplier'}) + $mod,$y1 + 20,($self->{'_flank'} * $self->{'_multiplier'}) + $key_box + $mod,$y1 + 20 + $key_box,$self->{'black'});
	$self->{'im'}->filledRectangle(($self->{'_flank'} * $self->{'_multiplier'}) + 1 + $mod,$y1 + 21,($self->{'_flank'} * $self->{'_multiplier'}) + $key_box - 1 + $mod,$y1 + 20 + $key_box - 1,$self->{'inside'});
	$self->write_text($self->{'_inside_label'},($self->{'_flank'} * $self->{'_multiplier'}) + $key_box + 5 + $mod,$y1 + 6 + $key_box);

	$mod += 150;

	$self->{'im'}->filledRectangle(($self->{'_flank'} * $self->{'_multiplier'}) + $mod,$y1 + 20,($self->{'_flank'} * $self->{'_multiplier'}) + $key_box + $mod,$y1 + 20 + $key_box,$self->{'black'});
	$self->{'im'}->filledRectangle(($self->{'_flank'} * $self->{'_multiplier'}) + 1 + $mod,$y1 + 21,($self->{'_flank'} * $self->{'_multiplier'}) + $key_box - 1 + $mod,$y1 + 20 + $key_box - 1,$self->{'outside'});
	$self->write_text($self->{'_outside_label'},($self->{'_flank'} * $self->{'_multiplier'}) + $key_box + 5 + $mod,$y1 + 6 + $key_box);

	$mod += 150;

	$self->{'im'}->filledRectangle(($self->{'_flank'} * $self->{'_multiplier'}) + $mod,$y1 + 20,($self->{'_flank'} * $self->{'_multiplier'}) + $key_box + $mod,$y1 + 20 + $key_box,$self->{'black'});
	$self->{'im'}->filledRectangle(($self->{'_flank'} * $self->{'_multiplier'}) + 1 + $mod,$y1 + 21,($self->{'_flank'} * $self->{'_multiplier'}) + $key_box - 1 + $mod,$y1 + 20 + $key_box - 1,$self->{'green'});
	$self->write_text($self->{'_unknown_label'},($self->{'_flank'} * $self->{'_multiplier'}) + $key_box + 5 + $mod,$y1 + 6 + $key_box);

	$mod += 150;

	$self->{'im'}->filledRectangle(($self->{'_flank'} * $self->{'_multiplier'}) + $mod,$y1 + 20,($self->{'_flank'} * $self->{'_multiplier'}) + $key_box + $mod,$y1 + 20 + $key_box,$self->{'black'});
	$self->{'im'}->filledRectangle(($self->{'_flank'} * $self->{'_multiplier'}) + 1 + $mod,$y1 + 21,($self->{'_flank'} * $self->{'_multiplier'}) + $key_box - 1 + $mod,$y1 + 20 + $key_box - 1,$self->{'dark_grey3'});
	$self->write_text($self->{'_helix_label'},($self->{'_flank'} * $self->{'_multiplier'}) + $key_box + 5 + $mod,$y1 + 6 + $key_box);
}

1;


=head1 NAME

LOCAL::lib::DrawMembrane.pm - draw a cartoon of an Alpha-helical transmembrane protein.

=head1 SYNOPSIS

use DrawTransmembraneSchematic;

 ## Create a hash and store topology, N-terminal and signal peptide information as follows

 my %topology;

 $topology{'TMHMM2.0'}{'topology'} = [5,22,46,68,89,111,135,157];
 $topology{'TMHMM2.0'}{'n_term'} = 'in';

 $topology{'PHOBIUS1.01'}{'signal'} = 22;
 $topology{'PHOBIUS1.01'}{'topology'} = [45,67,88,109,135,156];
 $topology{'PHOBIUS1.01'}{'n_term'} = 'out';

 $topology{'PRODIV0.92'}{'topology'} = [5,25,46,66];
 $topology{'PRODIV0.92'}{'n_term'} = 'out';

 $topology{'PROFPHD'}{'topology'} = [3,20,46,65,87,106,138,155,166,183];
 $topology{'PROFPHD'}{'n_term'} = 'in';

 $topology{'SIGNALP3.0'}{'signal'} = 24;
 $topology{'SIGNALP3.0'}{'topology'} = [];
 $topology{'SIGNALP3.0'}{'n_term'} = 'out';

 $topology{'MEMSAT3.0'}{'topology'} = [5,24,48,67,89,112,135,158];
 $topology{'MEMSAT3.0'}{'n_term'} = 'out';


 ## Title, sequence and the order to display topologies:

 my $title = 'DrawTransmembraneSchematic.pm Test';

 my $fa = 'MANMFALILVIATLVTGILWCVDGQSLNAPTSGNFFFAPVAPNPGWLVTGASVFPVLAIVLIVRSFIYEPFQIPSGSMMPTLLIGEVPMANMFFALILVIATLVTGILWCVDGQSLNAPTSGKFFFAPVKPKPGWLVTGASVFPVLAIVLIVRSFIYEPFQIPSGSMMPTLLIGDFILVEKFAYGIKDPIYQKTLIETGHPKRGDIVVFKYPEDPKLDYIKRAVGLPGDKVTYDPVSKELTIQPGCSSGQACENALPVTYSNVEPSDFVQTFSRRNGGEATSGFFEVPKNETKENGIRLSERKETLGDVTHRILTVPIAQDQVGMYYQQPGQQLATWIVPPGQYFMMGDNRDNSADSRYWGFVPEANLVGRATAIWMSFDKQEGEWPTGLRLSRIGGIH';

 my @order = ('MEMSAT3.0','PROFPHD','PRODIV0.92','PHOBIUS1.01','TMHMM2.0','SIGNALP3.0');


 ## If you want to use TTF fonts set this path and uncomment -ttf_font below

 my $font = '/usr/share/fonts/msttcorefonts/arial.ttf';

 my $im = DrawTransmembraneSchematic->new(-title=>$title,
                                          -order=>\@order,
                                         #-ttf_font=>$font,
                                          -sequence=>$fa,
                                          -topologies=>\%topology,
                                          -signal_overides=>1,
                                          -exclude_from_consensus=>'SIGNALP3.0');

 ## Or to add topologies as you go try this:

 my %hash;;
 $hash{'-method'} = 'MEMSAT4';
 $hash{'-topology'}= [5,24,48,67,89,112,135,158];
 $hash{'-n-term'} = 'out';
 $hash{'-signal'} = 0;

 $im->add_topology(%hash);


 ## Now write the image to a .png file

 open(OUTPUT, ">output.png");

 binmode OUTPUT;

 print OUTPUT $im->png;

 close OUTPUT;

=head1 DESCRIPTION

A module to draw an image showing multiple topology predictions for a protein and a consensus topology. It uses GD and allows the image to be written to a .png file.

The options are a set of tag/value pairs as follows:

  Option              Value                                         Default
  ------              -----                                         -------


   -topology_array         Array containing transmembrane helix          none
	                   boundaries

   -topology_string        Alternative to -topology, provide a string    none
                           containing the topology data in the form
		           A.11,31;B.41,59;C.86,107;D.145,166

   -n_term                 Location of the N-terminal of the sequence,   out
  	                   either 'in' or 'out'

   -signal                 Cleavage site for signal peptide              0

   -method                 Method name used with add_topology function   none

   -sequence               Protein sequence                              none

   -seq_length             Sequence length (not required if -sequence    500
                           is provided)

   -draw_kd                Draw Kyte-Doolittle plot                      1

   -kd_window              Kyte-Doolittle sliding window size            19

   -kd_multiplier          Height multiplier for Kyte-Doolittle plot     10

   -draw_key               Draw key                                      1

   -draw_grid              Draw grid lines every 100 residues            1

   -draw_consensus         Generate consensus topology                   1

   -title                  Title to add to the image                     none

   -inside_label           Label for inside loops                        Cytoplasmic

   -outside_label          Label for outside loops                       Extracellular

   -helix_label            Label for transmembrane helices               Transmembrane Helix

   -unknown_label          Label for unknown consensus regions           Unknown

   -signal_label           Label for signal peptide                      Signal Peptide,

   -kd_label               Label for Kyte-Doolittle plot                 Kyte-Doolittle

   -consensus_label        Label for consensus topology                  Consensus

   -topologies             Hash containing different topologies (see     none
                           above)

   -order                  Array containing different methods. This      none
                           will be the order the topologies are
                           displayed in

   -bar_height             Height of each track                          15

   -seperator              Gap between each track                        10

   -increment              Amount to increment after each topology       0

   -top                    Buffer at top of image                        40

   -multiplier             Factor to stretch image horizontally          1.6

   -text_buffer            Horizontal text buffer                        20

   -flank                  Buffer at either size of image                120

   -ttf_font               Path to TTF font, e.g.                        none
                           /usr/share/fonts/msttcorefonts/arial.ttf

   -ttf_font_size          Default size for TTF font. Use 7-9 with       8
                           Arial for best results

   -text_vertical_offset   Vertical text offset                          0

   -text_horizontal_offset Horizontal text offset                        0

   -signal_overides        Signal peptide overides consensus             1

   -exclude_from_consensus Don't include this method in the              none
                           consensus

   -dont_colour            Don't colour this method                      none

   -inside_rgb             Array containing RGB color code for          [255,255,255]
                           inside loops

   -outside_rgb            Array containing RGB color code for          [255,180,0]
                           outside loops

   -signal_rgb             Array containing RGB color code for          [243,48,130]
                           signal peptide


=head1 AUTHOR

Tim Nugent E<lt>timnugent@gmail.comE<gt>

=cut
