#!/usr/bin/perl

my @samples;
open FPIN, "<../resources/param.txt" or die;
while (<FPIN>) {
	if (/^sampleNames\t([^\n]+)\n$/) {
		my $line = $1;
		$line =~ s/^ //g;
		$line =~ s/ $//g;
		@samples = split / /, $line;
	}
}
close FPIN;


## build data structure

`mkdir ../output/exoncount/`;

my $ds;
foreach my $sample (@samples) {

	open FPIN, "<../output/".$sample."/exoncount/".$sample.".exoncount.markdup" or die;
	open FPOUT, ">../output/exoncount/".$sample.".exoncount.markdup" or die;
	while (<FPIN>) {
		print FPOUT $_;
	}	
	close FPIN;	
	close FPOUT;

	open FPIN, "<../output/".$sample."/exoncount/".$sample.".exoncount.rmdup" or die;
	open FPOUT, ">../output/exoncount/".$sample.".exoncount.rmdup" or die;
	while (<FPIN>) {
		print FPOUT $_;
	}	
	close FPIN;	
	close FPOUT;
}


