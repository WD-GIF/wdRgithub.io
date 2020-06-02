#!C:\usr\perl\perl.exe
use warnings;
use strict;
use 5.10.1;
use BioPerl;
use Bio::SeqIO;
use Bio::Seq;
use Data::Dumper;
my $alphabet;
my $id;
my $len;
my $seq;
my $in = Bio::SeqIO->new(-file => "final.fa", -format => 'Fasta', -alphabet =>"dna");
say "---------------------------------------------------第一题--------------------------------------";
say "序列依次对应下面primary_id的长度为:";
while (my$seqobj = $in->next_seq()){
	$id.=" ".$seqobj->id(); 	#id号
	 $seq.=" ".$seqobj->seq();	#序列
	#my $desc = $seqobj->desc();	#关于序列的描述信息
	$alphabet= $seqobj->alphabet();	#判断序列所属组学
	my $notelens=
	$len+= $seqobj->length();#每一段的序列长度

	print $seqobj->length(),"\t";
    }
say"";
 my $parnum=split/ /,$id;
 $parnum--;
 my $avelen=$len/$parnum;
 say "序列的primary_id分别为:",join "||",split/ /,$id;#以双竖线分割
 say "读取段数为:$parnum";
 say "平均长度为:$avelen";
##第二题,第一题产生的序列可用
say "---------------------------------------------------第二题--------------------------------------";
my @seqs=split " ",$seq;
#定义一个子程序来产生翻译序列
sub Trans{
	(my $dna)=shift;
 my $seqobj=Bio::Seq->new(-seq =>$dna, -alphabet =>'dna');
 return $seqobj->translate()->seq();
 }
 say"翻译为蛋白质序列 waiting……:";
say "Protein seq1".":".length(Trans($seqs[0])).":".Trans($seqs[0]);
say "Protein seq2".":".length(Trans($seqs[1])).":".Trans($seqs[1]);
say "Protein seq3".":".length(Trans($seqs[2])).":".Trans($seqs[2]);#输出也可以用Dumper
#第三题 统计CG含量,所提供的文件没有出现大小写的问题，都是大写，很方便
#第三题 同样利用第一题的$seq或者第二题的@seqs
say "---------------------------------------------------第三题--------------------------------------";
my $gcnum;
for(my $i=0;$i<scalar(@seqs)-1;$i++){
foreach (split"",$seqs[$i]){
	$gcnum++ if~/[GC]/;
}
}
$gcnum=sprintf "%4.2f",($gcnum/($avelen*3))*100;
say "GC含量约为:".$gcnum."%";
#将序列的互补序列翻译出来,依然利用地二三题的数据
say "---------------------------------------------------第四题--------------------------------------";
#先产生互补,定义一个小小的子程序
sub com{
	my $dna=shift;
	$dna=~tr/ACGTacgt/TGCAtgca/;
	$dna
}
#因为知道序列条数，所以我直接运行下面的代码，如果你不知道，请写一个循环来实现
my @revseqs=(com($seqs[0]),com($seqs[1]),com($seqs[2]));
say"互补序列翻译为蛋白质序列 waiting……:";#这次我用Dumper输出
my %protein=("蛋白1",Trans($revseqs[0]),"蛋白2",Trans($revseqs[1]),"蛋白3",Trans($revseqs[2]));
print Dumper(%protein);

########禁止转载抄袭,代码上传网址:https://github.com/WD-GIF/wdRgithub.io###############
#############----------王典的github--------------------------##########################





