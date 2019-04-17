use strict;
use warnings;
use Vcf;
use FaSlice;

my $fname = "path/to/input";
my $oname = "path/to/annotated/output";
my $fasta = "path/to/referece/fasa";
my $adj = 2; # how many adjacent sites to annotate singleton with?
my ($chr, $pos, $ref, $alt, $dp, $ns);

my %chr2fa = ();

sub getMotif {
  my $localseq=shift;
  my $adj=shift;
  my $subseq=$adj*2+1;

  my $altlocalseq = reverse $localseq;
  $altlocalseq  =~ tr/ACGT/TGCA/;

  my $ref1 = substr($localseq, $adj, 1);
  my $ref2 = substr($altlocalseq, $adj, 1);

  my $seqp;
  if($ref1 ~~ [qw( A C )]){
    $seqp = "$localseq\($altlocalseq\)";
  } else {
    $seqp = "$altlocalseq\($localseq\)";
  }

  return $seqp;
}

sub getType {
  my $ref=shift;
  my $alt=shift;
  my $adj=shift;
  my $seqp=shift;

  my $CAT = "${ref}${alt}";
  my $Category;
  if($CAT ~~ [qw( AC TG )]){ $Category = "AT_CG";}
  elsif($CAT ~~ [qw( AG TC )]){ $Category = "AT_GC";}
  elsif($CAT ~~ [qw( AT TA )]){ $Category = "AT_TA";}
  elsif($CAT ~~ [qw( GA CT )]){ $Category = "GC_AT";}
  elsif($CAT ~~ [qw( GC CG )]){ $Category = "GC_CG";}
  elsif($CAT ~~ [qw( GT CA )]){ $Category = "GC_TA";}

  if(substr($seqp, $adj, 2) eq "CG"){ $Category = "cpg_$Category";}
  return $Category;
}


open(my $fh, '<:encoding(UTF-8)', $fname)
  or die "Could not open file '$fname' $!";

open(my $ofh, '>', $oname) or die "Could not open file '$oname' $!";

print $ofh "POS\tREF\tALT\tMotif\tType\n";

while (<$fh>) {
    chomp;
    #"CHR", "POS", "REF","ALT","DP", "NS"
    ($chr, $pos, $ref, $alt, $dp, $ns) = split("\t");
    $alt = substr $alt, 0, 1;
  
    if (!exists($chr2fa{$chr})) {
        $chr2fa{$chr} = FaSlice->new(file=>$fasta, size=>100_000);
    }
    
    my $fa = $chr2fa{$chr};
    
    my $motif = $fa->get_slice($chr, $pos-$adj, $pos+$adj);
    $motif = getMotif($motif, $adj);
    
    my $type = getType($ref, $alt, $adj, $motif);
    
    print $ofh "$pos\t$ref\t$alt\t$motif\t$type\n";
}