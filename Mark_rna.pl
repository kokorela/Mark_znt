#!/usr/bin/perl

# command line example - perl mark_rna.pl chrN inputList

# 18    Cufflinks    exon    47390    48447    1000    -    .    gene_id "ENSG00000173213"; transcript_id "ENST00000308911"; exon_number "1"; FPKM "0.0176857771"; frac "1.000000"; conf_lo "0.000000"; conf_hi "0.060132"; cov "0.111308";
# ./init32_hsG chr18 78077248
# perl mark_rna.pl chr18 sample1/sample1_cuff/sample1_transcripts.gtf

$argc = scalar(@ARGV);
$chrCl = $ARGV[0];
$list = $ARGV[1];

chomp $list;

if($argc != 2) {
   print "argc $argc  chrCl $chrCl  list $list\n";
   exit;
}

$chrN = $chrCl . ".znt4";

unless(open(FINP,"+<$chrN")) {
   print "Unable to open chrN.znt4 file\n";
}
binmode FINP;
seek FINP,0,2;
$chrLen = tell(FINP)/4;
print "chrLen $chrLen\n";

unless(open(FLP,$list)) {
   print "Unable to open list file\n";
}

$fcp = 'mark_rna.intCode.' . $chrCl;
unless(open(FCP,">$fcp")) {
   print "Unable to open mark_rna.intCode.$chrCl $fcp file\n";
   exit;
}

$fep = 'mark_rnaErr.' . $chrCl;
unless(open(FEP,">$fep")) {
   print "Unable to open mark_rnaErr.$chrCl $fep file\n";
   exit;
}

$fp = 'mark_rnaOut.' . $chrCl;
unless(open(FP,"+>$fp")) {
   print "Unable to open mark_rnaOut.$chrCl $fp file\n";
   exit;
}
binmode FP;

$intCode = 1;

$smp1 = 0;
$smp2 = 0;
$smp3 = 0;
$smp4 = 0;
sub mkhs {
   for($n = $lpos; $n < $rpos; $n++) {
       $cp = $chrCl . '.' . $n;
       if(!$cph{$cp}) {
           if(!$isth{$tssIdVal}) {
               $isth{$tssIdVal} = $intCode;
               @tica = ($tssIdVal, $fpkmVal);
               $ich{$intCode} = [ @tica ];
               $cph{$cp} = $intCode;
               if($smp1 < 5) {
                   # print "11 cp $cp  cph_cp $cph{$cp}\n";
                   # print "11 intCode $intCode  @{ $ich{$intCode} }\n";
                   $smp1++;
               }
               $intCode++;
           }
           else {
               $tic = $isth{$tssIdVal};
               $cph{$cp} = $tic;
               if($smp2 < 5) {
                   # print "12 cp $cp  cph_cp $cph{$cp}\n";
                   # print "12 tic $tic  @{ $ich{$tic} }\n";
                   $smp2++;
               }
           }
       }
       else {
           if(!$isth{$tssIdVal}) {
               $isth{$tssIdVal} = $intCode;
               @tica = ($tssIdVal, $fpkmVal);
               $ich{$intCode} = [ @tica ];
               $cph{$cp} = $cph{$cp} . '|' . $intCode;
               if($smp3 < 5) {
                   # print "21 cp $cp  cph_cp $cph{$cp}\n";
                   # print "21 intCode $intCode  @{ $ich{$intCode} }\n";
                   $smp3++;
               }
               $intCode++;
           }
           else {
               $tic = $isth{$tssIdVal};
               $cph{$cp} = $cph{$cp} . '|' . $tic;
               if($smp4 < 5) {
                   # print "22 cp $cp  cph_cp $cph{$cp}\n";
                   # print "22 tic $tic  @{ $ich{$tic} }\n";
                   $smp4++;
               }
           }
       }
   }

   return;
}

$mp = 0;
while($temp = <FLP>) {
   chomp $temp;
   @ta = split (/\s/,$temp);
   $chrL = "chr" . $ta[0];
   $lpos = $ta[3];
   $rpos = $ta[4];
   $orient = $ta[6];
   $tssId = $ta[10];
   $tssIdTemp = $ta[11];
   $fpkm = $ta[14];
   $fpkmVal = $ta[15];
   if($chrCl eq $chrL) {
       if($temp =~ /exon/) {
           if($fpkm eq "FPKM" && $tssId eq "transcript_id") {
               $tssIdVal = $tssIdTemp;
               $tssIdVal =~ s/\"//;
               $tssIdVal =~ s/\"\;//;
               $fpkmVal =~ s/\"//;
               $fpkmVal =~ s/\"\;//;
               if($mp < 5) {
                   # print "chrL $chrL  lpos $lpos  rpos $rpos  tssIdTemp $tssIdTemp  tssIdVal $tssIdVal  fpkmVal $fpkmVal\n";
                   $mp++;
               }
               mkhs();
           }
           else {
               print "read phase error; fpkm $fpkm  transcript_id $tssId\n";
               exit;
           }
       }
   }
}

print "\n";
$mp1 = 0;
seek FP,0,0;
print FP pack "ic",0,0;
$fpSw = 1;
for $cp (keys %cph) {
   $ks = $cph{$cp};
   if(!$kslh{$ks}) {
       @ta = split(/\|/,$cph{$cp});
       $taLen = scalar(@ta);
       if($mp1 < 3 && $taLen > 2) {
           print "cp $cp  cph_cp $cph{$cp}  ks $ks\ntaLen $taLen ta @ta\n";
       }
       $fpEnd = $fpSw * 5;
       $kslh{$ks} = $fpEnd;
       for($n = 0; $n < $taLen - 1; $n++) {
           $chainChar = 1;
           print FP pack "ic",$ta[$n],$chainChar;
           print FEP "fpEnd $fpEnd  ta_n $ta[$n] $n   chainChar $chainChar\n";
           if($mp1 < 3 && $taLen > 2) {
               print "fpEnd $fpEnd  ta_n $ta[$n] $n   chainChar $chainChar\n";
           }
           $fpSw++;
       }
       $chainChar = 0;
       print FP pack "ic",$ta[$n],$chainChar;
       print FEP "fpEnd $fpEnd  ta_n $ta[$n] $n   chainChar $chainChar\n";
       if($mp1 < 3 && $taLen > 2) {
           print "fpEnd $fpEnd  ta_n $ta[$n] $n   chainChar $chainChar\n";
           $mp1++;
       }
      $fpSw++;
   }
}


print "\n";

$mp2 = 0;
seek FINP,0,0;
for($m = 0; $m < $chrLen; $m++) {
   $cp = $chrCl . '.' . $m;
   if($cph{$cp}) {
       $ks = $cph{$cp};
       $fpLoc = $kslh{$ks};
       if($mp2 < 5) {
           print "ks $ks  fpLoc $fpLoc\n";
           $mp2++;
       }
       seek FINP,4 * $m, 0;
       print FINP pack "i",$fpLoc;
   }
}

for $code (keys %ich) {
   print FCP "$code\t@{ $ich{$code} }\n";
}

$cphLen = scalar(keys %cph);
$ichLen = scalar(keys %ich);
$isthLen = scalar(keys %isth);
$kslhLen = scalar(keys %isth);

print "\nfpSw $fpSw  cphLen $cphLen  ichLen $ichLen  isthLen $isthLen  kslhLen $kslhLen\n";

close(FINP);
close(FLP);
close(FCP);
close(FEP);
close(FP);

exit;
