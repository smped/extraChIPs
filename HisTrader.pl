########################################################################################################
##                                                                                                    ##
##    HISTRADER: A tool to identify nucleosome free regions from ChIP-Seq of Histone Modifications    ##
##                                                                                                    ##
##                                                                                                    ##
##                                Written by Yifei Yan and Swneke D. Bailey                           ##
##                                   Copyright 2020 Swneke D. Bailey                                  ##
##                                                                                                    ##
########################################################################################################

use strict;
use warnings;
use POSIX;
use Getopt::Long;


# load genome fasta file into a hash / each sequence is one continuous string
# returns a hash
sub loadGenome{
   my $fname=shift;

   my $fh;
   my $name;
   my %h=();

   open($fh,"<",$fname) or die "Could not open $fname!\n";


   while(<$fh>){
      chomp();
      my $line=$_;

      if($line =~ /^\s*$/){
          next;
      }elsif($line =~ /^\s*#/){
          next;
      }elsif($line =~ /^>/){
          my @name=split(/\s+/, $line);
          $name=$name[0];
          $name=~ s/>//g; # remove fasta format
          next;
      }else{
          $h{$name} .= $line;
      }
   }

   close($fh);

   return(%h);
}

# extract fasta sequence given .bed file coordinates
sub getFasta{
    my $coord=shift; # a bed file of DHS coordinates -- change to a %hash
    my $genome_ref=shift; #hash containing genome sequence

    my @coord=split(/:|-/,$coord);
    my $seq;

    my $chr=$coord[0], my $start=$coord[1], my $end=$coord[2]; # redundant but clear

    if(exists ${$genome_ref}{$chr}){ # NOTE must use same naming convention
        my $pos=join("-",$start,$end);
        my $size=$end - $start;
        $seq=substr(${$genome_ref}{$chr},$start - 1,$size + 1); # CHECK 0 or 1 offset
    }

   return($seq);
}

sub movingAverageCentre{

  my $a = shift; # array with signal values;
  my $c = shift; # window fast;

  my $before =  floor((($c-1)/2));
  my $after =  ceil((($c-1)/2));

  my $tot = $#{$a} + 1;

  my @sum = @{$a}; # initialize sum
  my @count = (1) x $tot; #initialize count

  my $i = 1;

  while($i <= $before){
    my @new = @{$a};
    my @tmp = (0) x $i;
    splice(@new,$tot-$i,$i);
    splice(@new,0, 0,@tmp);

    for(my $j=0; $j<=$#{$a}; $j++){
       $sum[$j]+=$new[$j];
    }
    for(my $j=$i; $j<=$#{$a}; $j++){
       $count[$j]+=1;
    }
    $i+=1;
  }

  $i = 1;

  while($i <= $after){
    my @new = @{$a};
    my @tmp = (0) x $i;
    splice(@new,0,$i);
    splice(@new,$tot-$i, 0,@tmp);

    for(my $j=0; $j<=$#{$a}; $j++){
       $sum[$j]+=$new[$j];
    }
    for(my $j=0; $j<=$#{$a}-$i; $j++){
       $count[$j]+=1;
    }
    $i+=1;
  }

  my @movingAve=();
  for(my $j=0; $j<=$#{$a};$j++){
     push(@movingAve, $sum[$j]/$count[$j]);
  }

  return(@movingAve);
}

sub crossover{
  my $a=shift; # fast moving average
  my $b=shift; # slow moving average
  my $c=shift; # positions

  my @fast=@{$a};
  my @slow=@{$b};
  my @positions=@{$c};

  my @u=();
  for(my $i=0; $i<=$#fast; $i++){
     if($fast[$i] > $slow[$i]){
        push(@u,$i);
     }
  }

  my @pos=();
  for(my $i=0; $i<=$#u; $i++){
     push(@pos, $positions[$u[$i]]);
  }

  return(\@pos);
}

sub getNucFreeDiff{
    my $a = shift; #merge nucleosome

    my @dhs=();
    for(my $i=0; $i<=$#{$a}-1; $i++){
       my @start=split(/-/,${$a}[$i]);
       my @end=split(/-/,${$a}[$i+1]);
       my $tmp=join("-",$start[1],$end[0]);
       push(@dhs,$tmp);
     }

   return(@dhs);
}

sub round{
   my $a=shift;

   my $r=int($a + $a / abs($a * 2 || 1));

   return($r);
}

sub diff{
   my $a = shift;
   my $c = shift;

   #smooth signal
   my @ma = movingAverageCentre(\@{$a}, $c);
   my @diff=();
   for(my $i=0; $i<=$#{$a}-1; $i++){
      push(@diff, $ma[$i + 1] - $ma[$i]);
   }
   return(@diff);
}

sub diffOrd_2{
    my $a = shift;
    my $b = shift;

    my @o = ();
    for(my $i=0; $i <= $#{$a}; $i++){
        if(${$a}[$i] < 0){
           push(@o,$i + 2);
        }
    }
   my @positions=();
   for my $idx (@o){
       push(@positions, ${$b}[$idx]);
   }

   return(\@positions);
}

sub aveStep{
   my $a=shift;

   my $sum=0;
   my $n=$#{$a} + 1;
   for my $s (@{$a}){
     my @tmp=split(/-/,$s);
     my $size=$tmp[1] - $tmp[0];
     $sum+=$size;
   }
   my $average=($sum/$n);

   return($average);
}

sub my_min {
    my $min = shift;
    my $next;
    while(@_){
       $next = shift;
       $min = $next if($next < $min);
    }
    return $min;
}

sub increaseStep{
    my $a=shift; #signal
    my $b=shift; #step
    my $step=shift; #new step

    my @signal=@{$a};
    my @step=@{$b};
    my @dist=();

    for my $s (@step){
       my @s=split(/-/,$s);
       my $d=$s[1] - $s[0];
       push(@dist, $d);
    }

    my @newSig=();
    my @newStep=();

    my $new=0;
    my @tmpSig=();
    my @tmpStep=();
    for(my $i=0; $i<=$#dist; $i++){
      if($new < $step){ # will overshoot
         $new+=$dist[$i];
         my @pos=split(/-/,$step[$i]);
         push(@tmpSig, $signal[$i]);
         push(@tmpStep, @pos);
      }else{
         my $maxSig=my_max(@tmpSig);
         my $minStep=my_min(@tmpStep);
         my $maxStep=my_max(@tmpStep);
         my $region=join("-",$minStep,$maxStep);

         push(@newSig, $maxSig);
         push(@newStep, $region);

         @tmpSig=();
         @tmpStep=();

         $new=$dist[$i];
         push(@tmpSig, $signal[$i]);
         my @pos=split(/-/,$step[$i]);
         push(@tmpStep, @pos);

      }
    }
   return(\@newSig, \@newStep);
}

sub increaseSig{
    my $a=shift;
    my $step=shift;

    my @int=@{$a};

    my $end=$#int;
    my $e=0;
    my @new=();

    for(my $i=0; $i<=$end; $i+=$step){
       my @tmp=();
       if($e + $step <= $end){
         $e=$i+$step;
       }else{
         $e=$end;
       }
       for(my $j=$i; $j<$e; $j++){
         push(@tmp,$int[$j]);
       }
       my $max=my_max(@tmp);
       push(@new, $max);
   }
   return(\@new);
}


sub merge{
   my $a = shift;
   my $d = shift;

   # assumes sorted
   my @merge=();
   my @pos=split(/-/,${$a}[0]);
   my $start=$pos[0];
   my $end=$pos[1];

   my %merge=();

   $merge{$start}=$end;

   #print "@pos\n";

   for(my $i=1; $i <=$#{$a}; $i++){
      my @pos=split(/-/,${$a}[$i]);

      if($merge{$start} + $d >= $pos[0]){
         $merge{$start}=$pos[1];
      }else{
         $start=$pos[0];
         $end=$pos[1];
         $merge{$start}=$end;
      }
   }
   my @sorted=sortProbes(\%merge);
   for my $i (@sorted){
      my $nuc=join("-",$i, $merge{$i});
      push(@merge, $nuc);
   }
   return(@merge);
}

sub highRes{
    my $signal = shift;
    my $position = shift;

    my @hresSig=();
    my @hresPos=();


    for(my $i=0; $i<=$#{$position}; $i++){
       my @pos=split(/-/,${$position}[$i]);
       my $start=$pos[0];
       my $end=$pos[1];
       for(my $j=$start; $j<=$end - 1; $j++){
          my $pos=join("-",$j, $j+1);
          push(@hresSig, ${$signal}[$i]);
          push(@hresPos, $pos);
       }
    }
    return( \@hresSig, \@hresPos);
}

sub build_index{

   my $dataFile = shift;
   my $indexFile = shift;

   my $dh;
   my $ih;

   open($dh, "<", $dataFile) or die;
   open($ih, ">", $indexFile) or die;

   my $offset = 0;
   my $n=0;
   my %chr=();

   while (<$dh>){
       print $ih pack("Q>", $offset); #write index
       $offset = tell($dh);
       $n+=1;

       chomp();
       my @line=split(/\t/);

       if(!exists $chr{$line[0]}){ #obtain chromosome start line positions
           $chr{$line[0]}=$n;
       }
   }
   close($dh);
   close($ih);
   return(%chr);
}

# get a line from file using the index
sub line_with_index{
   my $dataFile = shift;
   my $indexFile = shift;
   my $lineNumber = shift;

   my $size;
   my $i_offset;
   my $entry;
   my $d_offset;

   $size = length(pack("Q>", 0));
   $i_offset = $size * ($lineNumber-1);

   seek($indexFile, $i_offset, 0) or return;
   read($indexFile, $entry, $size);

   $d_offset = unpack("Q>", $entry);

   seek($dataFile, $d_offset, 0);

   return scalar(<$dataFile>);
}


sub my_max{
    my $max = shift;
    my $next;
    while(@_){
        $next = shift;
        $max = $next if ($next > $max);
    }
    return $max;
}

sub getProbes{
    my $probes=shift;

    open(my $fh, "<", $probes) or die;

    my %probes=();

    while(<$fh>){
        chomp();

        my @line=split(/\t/);

        if($line[0] =~ /POSITION/ || $line[0] =~ /#/){
           next;
        }else{
           my $chr=$line[0];
           my $start=$line[1];
           my $probe=join("-",$line[1],$line[2]);
           if(!exists $probes{$chr}{$probe}){ # {chr}{start-stop};
              $probes{$chr}{$probe}=$start;
           }
        }
    }
    close($fh);

    return(%probes);
}

# sort using explicit subroutine name
sub sortProbes{
   my $probe_ref=shift;

   my @sorted = sort { ${$probe_ref}{$a} <=> ${$probe_ref}{$b} } keys %{$probe_ref};  # presuming numeric

   return( @sorted );
}


# need to keep track of position in array of probes and keep moving forward and not going back to the start!
sub getProbeInt{
    my $index=shift;
    my $ints=shift;
    my $bgChr_ref=shift;
    my $sortedProbes=shift;
    my $genome=shift;
    my $trim=shift;
    my $trimSize=shift;
    my $nucSize=shift;
    my $minSize=shift;
    my $t=shift;
    my $t2=shift;
    my $fixStep=shift;
    my $mergeDist=shift;
    my $method=shift;
    my $prefix=shift;
    my $outBG=shift;
    my $pThresh=shift;
    my $filter=shift;
    my %genome=();
    my $df;
    my $nf;
    my $bg;
    my $dFasta;
    my $nFasta;
    my $bgOut;
    my %sortedProbes=%$sortedProbes;
    my $pg=0;
    my $dhsBed=join(".",$prefix,"nfr.bed");
    my $nucBed=join(".",$prefix,"nuc.bed");
    my $missBed=join(".",$prefix,"missing.bed");

    #Print Output Filenames
    print "\nOutput Filenames:\nNucleosome Free Regions (NFRs) = $dhsBed\nNucleosome Occupied Regions = $nucBed\nNo NFRs Detected = $missBed\n";

    open(my $ih,"<",$index);
    open(my $dh,"<",$ints);

    open(my $db,">",$dhsBed);
    open(my $nb,">",$nucBed);
    open(my $mf,">",$missBed); # Sites without a NFR

   if($outBG == 1){
      $bgOut=join(".",$prefix,"bedGraph");;
      open($bg,">",$bgOut);
   }
   if(defined $genome){
       $pg=1;
       %genome=loadGenome($genome);
       $dFasta=join(".",$prefix,"nfr.fa");
       $nFasta=join(".",$prefix,"nuc.fa");
       print "DNA Sequences within NFRs = $dFasta\nDNA Sequences within occupied regions = $nFasta\n";
       open($df,">",$dFasta) or die;
       open($nf,">",$nFasta) or die;
    }

    my %int=();
    my %step=();

    my @s=sortProbes(\%{$bgChr_ref});

    for(my $i=0; $i <= $#s; $i++){

       if($s[$i] =~ /track/){ # If bedGraph has a track header
          next;
       }

       my $pos = ${$bgChr_ref}{$s[$i]};
       my @line=split(/\t/,line_with_index($dh, $ih, $pos));

       chomp(@line);

       my $chr=$s[$i];
       my $start=$line[1];
       my $end=$line[2];
       my $signal=$line[3];
       my $step=join("-",$start,$end);

       for my $probe ( @{$sortedProbes{$s[$i]}}){
           my @probe=split(/-/,$probe);

           if($probe[1] - $probe[0] >= $minSize){
              my @hresSig=();
              my @hresPos=();
              #while( $start < $probe[1]){
              while( $start <= $probe[1]){
                 @line=split(/\t/,line_with_index($dh, $ih, $pos));
                 chomp(@line);
                 $chr=$s[$i];
                 $start=$line[1];
                 $end=$line[2];
                 $signal=$line[3];
                 $step=join("-",$start,$end);
                 #check
                 if($chr =~ /[a-zA-Z]/){
                   if($chr eq $line[0]){
                     if($probe[1] >= $start && $probe[0] <= $end){
                      push(@{$int{$chr}{$probe}}, $signal);
                      push(@{$step{$chr}{$probe}}, $step);
                     }
                   }
                 }else{
                   if($chr == $line[0]){
                     if($probe[1] >= $start && $probe[0] <= $end){
                      push(@{$int{$chr}{$probe}}, $signal);
                      push(@{$step{$chr}{$probe}}, $step);
                    }
                   }
                 }
                 $pos+=1;
              }

              my $peak=join(":",$chr,$probe);

              if(@{$int{$chr}{$probe}} && @{$step{$chr}{$probe}}){

                 my ($highResSig, $highResStep)=highRes(\@{$int{$chr}{$probe}},\@{$step{$chr}{$probe}});
                 my ($newSignal, $newStep)=increaseStep(\@{$highResSig},\@{$highResStep},$fixStep);

                 @{$int{$chr}{$probe}}=();
                 @{$step{$chr}{$probe}}=();

                 if(@{$newSignal} && @{$newStep}){

                     my $maxSig=my_max(@{$newSignal});
                     my $cutoff=round($maxSig * $pThresh);

                     for(my $i=0; $i<=$#{$newSignal}; $i++){
                         if(${$newSignal}[$i] < $cutoff){
                            ${$newSignal}[$i] = 0;
                         }
                     }

                     my @dhs=();
                     my @merge=();
                     my @nuc=();
                     my @free=();
                     my @overlap=();
                     my @nucOver=();

                     if($method eq "DIFF" || $method eq "BOTH"){
                         my @first=diff(\@{$newSignal},$t);
                         my @second=diff(\@first, $t);
                         my $change=diffOrd_2(\@second,\@{$newStep});
                         my @pos=split(/-/,$probe);

                         if(@{$change}){
                           @merge=merge(\@{$change}, $mergeDist); #diff NUC
                           @dhs=getNucFreeDiff(\@merge); #diff DHS
                           if($method eq "DIFF"){
                              if(!@dhs){
                                 print $mf "$chr\t$pos[0]\t$pos[1]\n";
                              }
                           }
                         }else{
                           print $mf "$chr\t$pos[0]\t$pos[1]\tISSUE DIFF\n";
                         }

                     }
                     if($method eq "MA" || $method eq "BOTH"){
                         my @fast=movingAverageCentre(\@{$newSignal},$t);
                         my @slow=movingAverageCentre(\@{$newSignal},$t2);
                         my $move=crossover(\@fast,\@slow,\@{$newStep});
                         my @pos=split(/-/,$probe);
                         if(@{$move}){
                            @nuc=merge(\@{$move}, $mergeDist);
                            @free=getNucFreeDiff(\@nuc);
                            if($method eq "MA"){
                               if(!@free){
                                 print $mf "$chr\t$pos[0]\t$pos[1]\n";
                               }
                            }
                         }else{
                            print $mf "$chr\t$pos[0]\t$pos[1]\tISSUE MA\n";
                         }
                     }
                     if($method eq "BOTH"){
                         if(@dhs && @free){
                            @overlap=overlap_consensus(\@dhs,\@free);
                            @nucOver=overlap_union(\@merge, \@nuc);
                         }else{
                            my @pos=split(/-/,$probe);
                            print $mf "$chr\t$pos[0]\t$pos[1]\n";
                         }
                     }


                     if($method eq "DIFF"){
                         if(@dhs){
                            my @dhsF=filterBed(\@dhs,$filter);
                            printBed($chr,$probe,\@dhsF,$db); # print nucleosomes
                            printBed($chr,$probe,\@merge,$nb); # print nucleosomes
                            if(defined $genome){
                              printFasta($chr,$probe,\@dhsF,\%genome,$df,$trim,$trimSize);
                              printFasta($chr,$probe,\@merge,\%genome,$nf,$trim,$trimSize);
                            }
                         }
                     }

                     if($method eq "BOTH"){
                         if(@overlap){
                           my @overlapF=filterBed(\@overlap,$filter);
                           printBed($chr,$probe,\@overlapF,$db);
                           printBed($chr,$probe,\@nucOver,$nb);
                           if(defined $genome){
                              printFasta($chr,$probe,\@overlapF,\%genome,$df,$trim,$trimSize);
                              printFasta($chr,$probe,\@nucOver,\%genome,$nf,$trim,$trimSize);
                           }
                         }
                     }

                     if($method eq "MA"){
                        if(@free){
                          my @freeF=filterBed(\@free,$filter);
                          printBed($chr,$probe,\@freeF,$db);
                          printBed($chr,$probe,\@nuc,$nb);
                          if(defined $genome){
                             printFasta($chr,$probe,\@freeF,\%genome,$df,$trim,$trimSize);
                             printFasta($chr,$probe,\@nuc,\%genome,$nf,$trim,$trimSize);
                          }
                        }
                     }

                     # PRINT CONVERTED BEDGRAPH FILE
                     if($outBG == 1){
                        for(my $i=0; $i <= $#{$newSignal}; $i++){
                          my @pos=split(/-/,${$newStep}[$i]);
                          print $bg "$chr\t$pos[0]\t$pos[1]\t${$newSignal}[$i]\n";
                        }
                    }
                }
            }
        }
      }
    }
    close($ih);
    close($dh);
    close($db);
    close($nb);
    close($mf);
    # Close genome files if specified
    if(defined $genome){
       close($df);
       close($nf);
    }
    # Close bedgraph file if specified
    if($outBG==1){
      close($bg);
    }
}

sub filterBed{
   my $a=shift;
   my $b=shift; #filter size

   my @dInd=();

   for(my $i=$#{$a}; $i>=0; $i--){
       my @tmp=split(/-/,${$a}[$i]);
       my $size=$tmp[1] - $tmp[0];
       if($size > $b){
         push(@dInd, $i);
       }
   }
   for my $i (@dInd){
      splice(@{$a},$i,1);
   }
   return(@{$a});
}

sub printBed{
    my $a=shift; #chromosome
    my $b=shift; #peak
    my $c=shift; #array of coordinates
    my $d=shift; #file handle bedfile

    my $n=$#{$c} + 1;
    my $s=0;
    for my $i (@{$c}){
       $s+=1;
       my @tmp=split(/-/,$i);
         print $d "$a\t$tmp[0]\t$tmp[1]\t$a:$b\t$n\t$s\n"; #number of nfrs or nucleosomes
    }
}

sub printFasta{
    my $a=shift; #chromosome
    my $b=shift; #peak
    my $c=shift; #array of coordinates
    my $d=shift; #genome
    my $e=shift; #file handle fasta
    my $f=shift; #trim;
    my $g=shift; #trimSize

    my $n=$#{$c} + 1;
    my $s=0;

    for my $i (@{$c}){
       my $coord;
       $s+=1;
       my @tmp=split(/-/, $i);

         if($f == 1){
           my $mid=round(($tmp[1] + $tmp[0])/2);
           my $ext=round($g/2);
           my $start=$mid - $ext;
           my $end=$mid + $ext;
           my $region=join("-",$start,$end);
           $coord=join(":",$a,$region);
         }else{
           $coord=join(":",$a,$i);
         }

         my $seq=getFasta($coord,\%{$d});
         print $e ">$coord\_$b\_$n\_$s\n$seq\n";
    }
}

sub printHeader{
    print "\n\n";
    print "########################################################################################################\n";
    print "##                                                                                                    ##\n";
    print "##    HISTRADER: A tool to identify nucleosome free regions from ChIP-Seq of Histone Modifications    ##\n";
    print "##                                                                                                    ##\n";
    print "##                                                                                                    ##\n";
    print "##                                Written by Yifei Yan and Swneke D. Bailey                           ##\n";
    print "##                                   Copyright 2020 Swneke D. Bailey                                  ##\n";
    print "##                                                                                                    ##\n";
    print "########################################################################################################\n";
}

sub overlap_union{
    my $a=shift;
    my $b=shift;

    my %regions=();
    my %overlap=();
    my @overlap=();
    my %dups=();

    for my $s (@{$a}, @{$b}){
        my @t=split(/-/,$s);
        $regions{$s}=$t[0]; #removes duplicates!!!
        $dups{$s}+=1;
    }

    my @regions=sortProbes(\%regions);
    #my $np=$regions[0];

    for(my $i=0; $i<=$#regions; $i++){
        for(my $j=0; $j<=$#regions; $j++){

           if($i != $j){
             #my @s=split(/-/,$np);
             my @s=split(/-/,$regions[$i]);
             my @t=split(/-/,$regions[$j]);

              if($s[1] >= $t[0] && $s[0] <= $t[1]){
                my @tmp=($s[0],$t[0]);
             #my $start=$t[0];
                my $start=my_min(@tmp);
                my @tmp2=($s[1],$t[1]);
                my $end=my_max(@tmp2);
                my $np=join("-",$start,$end);
                 $overlap{$start}=$np;
              }#else{
              # my $np=$regions[$i];
              # $overlap{$s[0]}=$np;
              #}
           }
        }
    }

    for my $d (sort keys %dups){
       if($dups{$d} > 1){
         push(@overlap, $d);
       }
    }

    for my $o (sort keys %overlap){
       push(@overlap, $overlap{$o});
    }

    %regions=();
    for my $r (@overlap){
       my @tmp=split(/-/, $r);
       $regions{$r}=$tmp[0];
    }
    @overlap=();
    @overlap=sortProbes(\%regions);
    return(@overlap);
}

sub overlap_consensus{
    my $a=shift;
    my $b=shift;

    my %regions=();
    my %overlap=();
    my @overlap=();
    my %dups=();

    for my $s (@{$a}, @{$b}){
        my @t=split(/-/,$s);
        $regions{$s}=$t[0]; #removes duplicates!!!
        $dups{$s}+=1;
    }

    my @regions=sortProbes(\%regions);
    #my $np=$regions[0];

    for(my $i=0; $i<=$#regions; $i++){
        for(my $j=0; $j<=$#regions; $j++){

           if($i != $j){
           #my @s=split(/-/,$np);
           my @s=split(/-/,$regions[$i]);
           my @t=split(/-/,$regions[$j]);

           #if($s[1] >= $t[0] && $s[0] <= $t[1]){
           if($s[1] > $t[0] && $s[0] < $t[1]){
             my @tmp=($s[0],$t[0]);
             #my $start=$t[0];
             my $start=my_max(@tmp);
             my @tmp2=($s[1],$t[1]);
             my $end=my_min(@tmp2);
             my $np=join("-",$start,$end);
             $overlap{$start}=$np;
           }#else{
         #    $np=$regions[$i];
         #  }
           }
        }
    }

    for my $d (sort keys %dups){
       if($dups{$d} > 1){
         push(@overlap, $d);
       }
    }

    for my $o (sort keys %overlap){
       push(@overlap, $overlap{$o});
    }

    %regions=();
    for my $r (@overlap){
       my @tmp=split(/-/, $r);
       $regions{$r}=$tmp[0];
    }
    @overlap=();
    @overlap=sortProbes(\%regions);
    return(@overlap);
}

### MAIN PROGRAM
### current the main program is getProbeInt ...
### the command line flags/options

my $chip;
my $probes;
my $genome;
my $method="BOTH";
my $fixStep=25;
my $mergeMulti=3; # multiplier for mergeing 3*25 = 75bp or 1/2 nucleosome - moving average
my $maMulti=3;
my $minPeak=500;
my $nucSize=150;
my $trim=0;
my $trimSize=100;
my $help=0;
my $pMax=0;
my $filter=1000;
my $outBG=0;
my $outPrefix="Histrader";


GetOptions( "bedGraph:s"   => \$chip,   # --bedGraph ... Histone ChIP-seq bedgraph
	    "peaks:s"      => \$probes,   # --peaks ...  Bed file
	    "genome:s"     => \$genome,   # --genome ... Genome fasta (Optional)
            "method:s"	   => \$method,   # --method ... method to call valleys (Optional)
            "step:i"       => \$fixStep,  # --step ... fixed bedgraph step (Optional)
            "out:s"        => \$outPrefix, # outputFile prefix
            "minSize:i"    => \$minPeak, # --minSize ... minimum peak size (Optional)
            "mergeMulti:i" => \$mergeMulti, # multiplier for merging
            "maMulti:i"    => \$mergeMulti, # multiplier for merging
            "nucSize:i"    => \$nucSize, # nucleosome size in bp
            "trim"	   => \$trim,
            "trimSize:i"   => \$trimSize, # size in bp for output fasta file
            "outBG"        => \$outBG,
            "pMax:f"         => \$pMax, #percentage of max in peak
            "filter:i"       => \$filter, #percentage of max in peak
            "help"         => \$help); #--help ... Print help message

printHeader;

if(defined $chip && defined $probes){
  print "\n\nIdentifying Valleys in $chip at positions specified in $probes\n\n";
  print "Parameters Used:\nBedGraph File = $chip\nPeak File = $probes\nOutput File Prefix = $outPrefix\n";
  if(defined $genome){
    print "Genome File = $genome\n";
  }else{
    print "Genome File = NOT SPECIFIED\n";
  }
  print "Method Used = $method\nStep Used = $fixStep\nMinimum Peak Size = $minPeak\nNucleosome Size = $nucSize\n";
  if($trim == 1){
     print "Trimming = YES\n";
     print "Trim Size = $trimSize\n";
  }else{
     print "Trimming = NO\n";
     print "Trim Size = NOT SPECIFIED\n";
  }
  if($outBG == 1){
     print "Output BedGraph = YES\n";
  }else{
     print "Output BedGraph = NO\n";
  }
}else{
  $help=1;
}

if($help==1){
   print "\nUSAGE: perl HISTRADER.pl --bedGraph ChIP.bedGraph --peaks ChIP.bed\n";
   print "\nOPTIONS:\n\n";
   print "--bedGraph\tSpecify ChIP-Seq signal file (bedGraph format) (required)\n\n";
   print "--peaks\t\tSpecify the BROAD PEAK file (bed format) (required)\n\n";
   print "--genome\tSpecify genome fasta file (optional)\n\t\tUsed to extract the DNA sequence within the valleys/footprints/NFRs\n\n";
   print "--trim\t\tSpecify that fasta sequences should be trimmed\n\t\tDefault (OFF).\tREQUIRES --genome and --trimSize.\n\n";
   print "--trimSize\tThe length of fasta sequence in base pairs (bp) returned for each valley / footprint / NFR. Sequences are trimmed from the centre of each valley / footprint / NFR.\n";
   print "\t\tREQUIRES --trim and --genome\n\n";
   print "--out\t\tThe output file prefix\n\t\tDefault=Histrader (optional)\n\n";
   print "--method\tMethod used to call valleys / footprints/ NFRs\n\t\tOptions = MA (moving average) / DIFF (differencing) / BOTH (merge of MA and DIFF)\n\t\tDefault = BOTH\n\n";
   print "--step\t\tThe fixed step size to use. The ChIP-seq profile will be converted to have fixed step equal to this number in base pairs (bp).\n\t\tDefault=25 (optional)\n\n";
   print "--minSize\tThe mininum peak size. Valleys/Footprints/NFRs will NOT be called in peaks that are smaller than this value\n\t\tDefault=500 (optional)\n\n";
   print "--nucSize\tThe estimated nucleosome size in base pairs (bp).This value is used for the moving average and smoothing the signal.\n\t\tThis value should be divisible by the step size specified (--step).\n\t\tDefault: 150 (~147).\t150 / 25 = 6 (6 bins is equivalent to 1 nucleosome)\n\n";
   print "--mergeMulti\tThe multiplier of step (--step) to use for merging.\n\t\tDefault = 3\t( 3 x 25 = 75 or ~ 1/2 a nucleosome)\n\n";
   print "--maMulti\tThe moving average multiplier for the slow moving average.\n\t\tDefault = 3 (ie. Slow moving average is equilavent to 3 nucleosomes [3 x 150])\n\n";
   print "--pMax\t\tFraction of max peak height to use as a threshold. Use to reduce calls in regions with low signal.\n\t\tFor example, use 0.25 (25%) to exclude signal less than 25% of the current peak's max height.\n\t\tDefault = 0 (or 0%) (Optional)\n\n";
   print "--filter\tFilter NFRs greater than this value. Needed for --pMAX, which can lead to large regions called as NFRs depending on the inputted peaks\n\t\tDefault = 1000 (bp) (Optional).\n\n";
   print "--outBG\t\tOutput the fixed step ChIP-Seq signal at the specified broad peaks (bedGraph format).\n\n";
   print "--help\t\t(prints this message)\n\n";
   exit();
}


my $mergeDist=$fixStep*$mergeMulti;
my $tmp=$nucSize/$fixStep; #number of bins per nucleosome
my $t=round($tmp);
$tmp=($nucSize*$maMulti)/$fixStep;
my $t2=round($tmp);

# print Additional parameters
print "Merge Distance = $mergeDist\nBins Per Nucleosome = $t\n";
print "pMax = $pMax\nFilter = $filter\n";
my %probes=getProbes($probes);

my %sortedProbes=();
for my $chr (keys %probes){
   my @s=sortProbes(\%{$probes{$chr}});
   $sortedProbes{$chr}=\@s;
}
%probes=();

my $index=join(".",$outPrefix,"idx");
my %bgChr=build_index($chip,$index);
getProbeInt($index,$chip,\%bgChr,\%sortedProbes,$genome,$trim,$trimSize,$nucSize, $minPeak, $t, $t2, $fixStep, $mergeDist, $method, $outPrefix, $outBG, $pMax, $filter);

print "\nFINISHED!\n";
