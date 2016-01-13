(* feature.ml *)

open Common
open Collection

type strand = | Forward  
              | Reverse
	      | Unstranded
	      | Unknown


let strand_of_string s =
  match s with
    | "++"
    | "+" -> Forward
    | "+-"
    | "-" -> Reverse
    | "." -> Unstranded
    | "?" -> Unknown
    |  _  -> failwith "strand_of_string : the strand should be provided as the 7th field of your gtf file\n" 


let strand_to_string = function
  | Forward -> "+"
  | Reverse -> "-"
  | Unstranded -> "."
  | Unknown -> "?"


let strand_to_string2 = function
  | Forward -> "p"
  | Reverse -> "m"
  | Unstranded -> "."
  | Unknown -> "?"



(* Takes a string v and checks whether it represents a fdd value field = be an integer, is greater than 10 and even *)
let is_gff_value_field v =
  let i = int_of_string v in
    ((i>=10) && ((i mod 2)==0));;
      
  
(* From a string of what we hope are integers that are more than 10 and even separated by commas,
   get a list of such integers.
   Raise a failure if one of those conditions is not met since these values would not then 
   correspond to gff value fields. *)
let get_list_of_gff_value_field s =
  let lval = split ',' s in
    if (List.for_all is_gff_value_field lval) then
      lval
    else
      begin
	failwith "with the -m n option and when n is not an integer, n has to be a list of integers greater than 10 and even, separated by commas\n"; 
	[]
      end;;
    
    
let to_non_null_string v=
    if (v="") then
      "."
    else
      v

(* checks whether the value v of a gff (key,value) pair is bordered by two double quotes and a semicolon 
   and if this is not the case adds some *)
let value_to_ucsc_value v = 
  let n=String.length v in
    try
      (
	if (((String.get v 0)='"') && ((String.sub v (n-2) 2)="\";")) then
	  v
	else
	  ("\"")^(v)^("\";")
      )
    with 
	Invalid_argument _ -> ("\"")^(v)^("\";")

(* Corresponding to a personalized type of gff format with tabulations first and then a list of attributes
   formatted as (key,value).*)
module Feature =
struct
  type t = {
    seq: string;                      (* in 1st field in gffflex *)
    start: int;                       (* in 4th field in gffflex *)
    stop: int;                        (* in 5th field in gffflex *)
    str: strand;                      (* in 7th field in gffflex *)
    ftype: string;                    (* in 3rd field in gffflex *)
    source: string;                   (* in 2nd field in gffflex *)
    score: string;                    (* in 6th field in gffflex *)
    frame: string;                    (* in 8th field in gffflex
					 in fact could only be . or 0 or 1 or 2, and these 3 only in case of CDS feature 
				      *)
    attlist: (string * string) list   (* in 9th field in gffflex
					 a list of attributes of the form (key, value) 
				      *) 
  }

  let create seq sta sto str ft src sc fr al = 
    {seq=seq; start=sta; stop=sto; str=str; ftype=ft; source=src; score=sc; frame=fr; attlist=al}
  
  let null = create "" (-1) (-1) Unstranded "null" "null" "." "." []
  let isnull f = (f=null)
    
  let compare_coord f1 f2 = 
    if ((Pervasives.compare f1.seq f2.seq)!=0) then
      Pervasives.compare f1.seq f2.seq
    else
      begin
	if ((Pervasives.compare f1.start f2.start)!=0) then
	  Pervasives.compare f1.start f2.start
	else
	  Pervasives.compare f1.stop f2.stop
      end

  let interval f = (f.start,f.stop)   
  let setattlist al f = {f with attlist=al}

  (* makes a string for gff flex output from a feature attribute list
     note that it takes as input the ucsc boolean option, which if it is true
     is the fact of adding two double quotes and one semi colon to each value that does not already have them 
  *)
  let string_attlist_in_gff_flex ucsc f = 
    let s= ref "" in
      (* note that here we are adding a space at the end of each line, that is why we do clean_end_string *)
      if (ucsc) then
	List.iter (fun (k,v) -> s:=(!s)^(" ")^(k)^(" ")^(value_to_ucsc_value (to_non_null_string v))) f.attlist
      else
	List.iter (fun (k,v) -> s:=(!s)^(" ")^(k)^(" ")^((to_non_null_string v))) f.attlist;
      (clean_end_string !s)
	
	    

  (* print_gff_flex takes as input:
     - o which is the output channel corresponding to the file where we need to write file1 with overlap info 
     - u which is the boolean ucsc printing option = we add two double quotes and one semi colon to all values of a pair (key,value)
       when this is not already the case (will modify all the values of the output file, even those that were not generated by overlap)
     - f which is a feature
     and prints in o the feature in the gff flexible format. *)
  let print_gff_flex o ucsc f =
    Printf.fprintf o "%s\t%s\t%s\t%i\t%i\t%s\t%s\t%s\t%s\n" f.seq f.source f.ftype f.start f.stop f.score (strand_to_string f.str) f.frame (string_attlist_in_gff_flex ucsc f)

  let isplus f = match f.str with
    | Forward -> true
    | _ -> false
  let isminus f = match f.str with
    | Reverse -> true
    | _ -> false
  let isplusorminus f = match f.str with
    | Forward | Reverse -> true
    | _ -> false
  let issamestr str f = (f.str=str)
    (* having different strands is only valid for + or - strand features *)
  let isdiffstr str f = (((str=Forward)||(str=Reverse))&&(isplusorminus f)&&(f.str!=str))
  
  let interval f = (f.start,f.stop)

  let intersection f1 f2 =
    if (f1.seq=f2.seq) then
      ((max f1.start f2.start), (min f1.stop f2.stop))
    else
      raise (Invalid_argument "intersection")

  (* f1 generally (not stricly) included in f2 *)
  let inclusion f1 f2 =
    if (f1.seq=f2.seq) then
      Common.finclusion (interval f1) (interval f2)
    else
      raise (Invalid_argument "inclusion")
  
    
  let seq f = f.seq
  let start f = f.start
  let stop f = f.stop
  let str f = f.str
  let ftype f = f.ftype
  let source f = f.source
  let score f = f.score
  let frame f = f.frame
  let attlist f = f.attlist
end



module Exon =
struct
  type t = {
	chrom: string;
	gbeg: int;
	gend: int;
	str: strand;
	score: float;  
 	cat: string;   (* 2nd field in gff (used to be 3rd field for makeSP but does not make sense since always exon then *)
	trid: string;
	gnid: string;
  }

  let create chr gb ge st sc ct tr gn = { chrom=chr; gbeg=gb; gend=ge; str=st; score= sc; cat=ct; trid=tr; gnid=gn}
  let null = create "" (-1) (-1) Forward  0.0 "" "" ""
  let isnull e = (e=null)

  (* when we compare two exons we compare them by gene first and then by coord, independantly of the strand 
     was here for makesp but many other ways to compare, some being listed below *)
  let compare e1 e2 = 
    if ((Pervasives.compare e1.gnid e2.gnid) !=0) then
      Pervasives.compare e1.gnid e2.gnid
    else
      begin
	if ((Pervasives.compare e1.gbeg e2.gbeg) !=0) then
	  Pervasives.compare e1.gbeg e2.gbeg
	else
	  Pervasives.compare e1.gend e2.gend
      end
  
  (* assuming that two exons are on the same strand (and that strand is + or -), compare their two chr
     and in case they are different return the alphabetical order for those. In case the chr are identical, 
     will return 1 when e1 is on 5' of e2, -1 when e1 is on 3' of e2 and 0 if they are equally on 5'. *)
  let compare_pos e1 e2 =
    if ((Pervasives.compare (e1.chrom) (e2.chrom))!=0) then
      Pervasives.compare e1.chrom e2.chrom
    else
      begin
	match e1.str with
	  | Reverse -> Pervasives.compare e2.gend e1.gend
	  | _ -> Pervasives.compare e1.gbeg e2.gbeg   (* might not be the ideal thing to do for not forward strand, to refine later on if needed *)
      end

  (* assuming that two exons are on the same strand (and that strand is + or -), compare their two chr
     and in case they are different return the alphabetical order for those. In case the chr are identical, 
     will return 1 when e1 is genomically before e2, -1 when e1 is genomically after e2 and 0 if they are 
     identical genomically. This is done independantly of the strand. *)
  let compare_pos_gx e1 e2 =
    if ((Pervasives.compare (e1.chrom) (e2.chrom))!=0) then
      Pervasives.compare e1.chrom e2.chrom
    else
      begin
	if ((Pervasives.compare e1.gbeg e2.gbeg) !=0) then
	  Pervasives.compare e1.gbeg e2.gbeg
	else
	  Pervasives.compare e1.gend e2.gend
      end


  (* given an output channel, a ucsc boolean and an exon, prints the exon in gtf format *)
  let print_gtf o u e = 
    if (isnull e) then
      ()
    else
      begin
	if u then
	  Printf.fprintf o "%s\t%s\texon\t%i\t%i\t%f\t%s\t.\tgene_id %s transcript_id %s\n" e.chrom e.cat e.gbeg e.gend e.score (strand_to_string e.str) (value_to_ucsc_value e.gnid) (value_to_ucsc_value e.trid)
	else
	  Printf.fprintf o "%s\t%s\texon\t%i\t%i\t%f\t%s\t.\tgene_id %s transcript_id %s\n" e.chrom e.cat e.gbeg e.gend e.score (strand_to_string e.str) e.gnid e.trid
      end	  


  let chrom e = e.chrom
  let gbeg e = e.gbeg
  let gend e = e.gend
  let str e = e.str
  let score e = e.score
  let cat e = e.cat
  let trid e = e.trid
  let gnid e = e.gnid

  let settrgn tr gn e = {e with trid=tr; gnid=gn}
end 


(* this function takes as input a feature supposed to be an exon extracted from gff with a gene in 10th field
   and a transcript in 12th field, and makes an exon object out of it. It will raise an exception if the 11th
   field was not transcript_id, and this error will be catched by the function Input.read_gff_into_tr_list
*)
let feature_to_exon f =
  if ((fst (List.nth (Feature.attlist f) 1)) = "transcript_id") then 
    Exon.create 
      (Feature.seq f) 
      (Feature.start f) 
      (Feature.stop f) 
      (Feature.str f) 
      (try (float_of_string (Feature.score f)) with | Failure _ -> 0.0)
      (Feature.source f) 
      (snd (List.nth (Feature.attlist f) 1)) 
      (snd (List.nth (Feature.attlist f) 0))
  else
    raise Not_found


(* this function takes as input an array of exons that have been put together by a projection 
   and makes a single exon out of them.
   I should catch an exception for exarr.(0), not done here.*)
let exarr_to_projex exarr =
  let e0= exarr.(0) and lex = Array.to_list exarr in
    Exon.create 
      (Exon.chrom e0)
      (List.fold_left min (Exon.gbeg e0) (List.map Exon.gbeg lex))
      (List.fold_left max (Exon.gend e0) (List.map Exon.gend lex))
      (Exon.str e0)
      (Common.avg_float (List.map Exon.score lex))
      (Exon.cat e0)
      "."
      "."
      


module Transcript =
struct
  type t = {
	chrom: string;
	gbeg: int;   
	gend: int;   (* gbeg is always smaller than gend, since those are genomic beg and end, not 5' and 3' extremities *)
	str: strand;
	cat: string;
	exlist: Exon.t list;
	trid: string;
	gnid: string;
  }
	  
  let create chr gb ge st ct el tr gn = 
    { chrom=chr; gbeg=gb; gend=ge; str=st; cat=ct;  exlist=el; trid=tr; gnid=gn}
  let null = create "" (-1) (-1) Forward "" [] "" ""
  let isnull t = (t=null)

  (* when we compare two transcripts we compare them by gene first and then by coord, independantly of the strand 
     was here for makesp for exons so needed to follow here, but many other ways to compare, some being listed below *)
  let compare t1 t2 = 
    if ((Pervasives.compare t1.gnid t2.gnid) !=0) then
      Pervasives.compare t1.gnid t2.gnid
    else
      begin
	if ((Pervasives.compare t1.gbeg t2.gbeg) !=0) then
	  Pervasives.compare t1.gbeg t2.gbeg
	else
	  Pervasives.compare t1.gend t2.gend
      end
  
  (* compare two transcripts according to strand and then to 5' to 3' position (without assuming they are on the same chr, but taking this info into account). 
     If they are on different strands then the lower one is the one on the + and if they are on the same strand then if this strand is +
     we compare their beg and returns the result of this comparison, otherwise we compare their end and take the opposite of
     the result of this comparison. *)
  let compare_strand_and_pos t1 t2 =
    match (t1.str,t2.str) with
      | (Forward,Reverse) -> -1
      | (Reverse,Forward) -> 1
      | (Forward,Forward) -> if ((Pervasives.compare (t1.chrom) (t2.chrom))!=0) then Pervasives.compare t1.chrom t2.chrom else Pervasives.compare t1.gbeg t2.gbeg
      | (Reverse,Reverse) -> if ((Pervasives.compare (t1.chrom) (t2.chrom))!=0) then Pervasives.compare t1.chrom t2.chrom else Pervasives.compare t2.gend t1.gend
      | (Forward,_) ->  -1
      | (_,Reverse) -> 1
      | _ -> if ((Pervasives.compare (t1.chrom) (t2.chrom))!=0) then Pervasives.compare t1.chrom t2.chrom else Pervasives.compare t1.gbeg t2.gbeg  (* might not be the ideal thing to do, to refine later on if needed *)

  (* compare two transcripts according to strand and genomic position (without assuming they are on the same chr, but taking this info into account). 
     note: for tr on the same strand and chromosome, the result does not depend on the strand since order is gx and not from 5' to 3' *)
  let compare_strand_and_pos_gx t1 t2 =
    match (t1.str,t2.str) with
      | (Forward,Reverse) -> -1
      | (Reverse,Forward) -> 1
      | (Forward,Forward) | (Reverse,Reverse) -> if ((Pervasives.compare (t1.chrom) (t2.chrom))!=0) then Pervasives.compare t1.chrom t2.chrom else Pervasives.compare t1.gbeg t2.gbeg
      | (Forward,_) ->  -1
      | (_,Reverse) -> 1
      | _ -> if ((Pervasives.compare (t1.chrom) (t2.chrom))!=0) then Pervasives.compare t1.chrom t2.chrom else Pervasives.compare t1.gbeg t2.gbeg  (* might not be the ideal thing to do, to refine later on if needed *)

  (* reports whether a transcript is monoexonic or not *)
  let monoex t = ((List.length t.exlist)==1)

  (* reports whether a transcript is spliced or not *)
  let splicedtr t = ((List.length t.exlist)>=2)

  (* reports whether two transcripts have overlapping coordinates (do not look at strand) *)
  let overlap t1 t2 =
    (((Pervasives.compare (t1.chrom) (t2.chrom))==0) && (Common.foverlap (t1.gbeg,t1.gend) (t2.gbeg,t2.gend)))

  let same_strand_and_chrom t1 t2 =
     (((t1.str)==(t2.str)) && ((Pervasives.compare (t1.chrom) (t2.chrom))==0))

  let interval t = (t.gbeg,t.gend)   

  (* reports the intron list of a transcript in the 5' to 3' order, not assuming that the exon list
     of the tr is itself ordered from 5' to 3' (means sort it beforehand). Introns are simple intervals. *)
  let intronlist t =
    let lex = List.sort Exon.compare_pos t.exlist in
    let rec auxnorm le =
      match le with
	| e1::(e2::q as q2) -> (((Exon.gend e1)+1),((Exon.gbeg e2)-1))::auxnorm(q2)
	| _ -> []
    in
    let	rec auxback le =
      match le with
	| e1::(e2::q as q2) -> (((Exon.gend e2)+1),((Exon.gbeg e1)-1))::auxback(q2)
	| _ -> []
    in
      match t.str with
	| Reverse -> auxback lex
	| _ -> auxnorm lex;;


  (* takes two spliced transcripts which are overlapping as input (not checked but would be better), 
     and reports whether overlapping introns of those are compatible (here means the same but in 
     the future we might want to allow for some nt of differences). Here I am using the most naive method
     which consists in making all the possible comparisons between t1 and t2 introns (symmetrical). *)
  let compatible_introns t1 t2 =
    let il1 = intronlist t1 and il2 = intronlist t2 in
    let lintrpair = ref [] in
    let u = List.iter (fun i1 -> List.iter (fun i2 -> (lintrpair:=(i1,i2)::(!lintrpair))) il2) il1 in
    let loverintrpair = List.filter (fun (i1,i2) -> if (Common.foverlap i1 i2) then true else false) (!lintrpair) in
    List.fold_left (fun bool (i1,i2) -> (bool && (i1=i2))) true loverintrpair;;

  (* this was added on 03/29/2012 following a sugegstion by Ben Brown: in order not to merge transcripts 
     that have the 3' utr of one overlapping the 5' utr of the second one, I just need to check that
     the two transcripts have at least one intron in common. This is the above function but with
     additional check that there is a pair of overlapping introns. The above function might disappear in future? *)
  let have_overlapping_and_compatible_introns t1 t2 =
    let il1 = intronlist t1 and il2 = intronlist t2 in
    let lintrpair = ref [] in
    let u = List.iter (fun i1 -> List.iter (fun i2 -> (lintrpair:=(i1,i2)::(!lintrpair))) il2) il1 in
    let loverintrpair = List.filter (fun (i1,i2) -> if (Common.foverlap i1 i2) then true else false) (!lintrpair) in
    (loverintrpair!=[]) && (List.fold_left (fun bool (i1,i2) -> (bool && (i1=i2))) true loverintrpair);;

  (* requires two spliced transcripts which are overlapping as input (not checked but would be better), 
     and reports whether there is no exon of one that overlaps even by 1 bp an intron of the other. 
     Here I am using the most naive method which consists in making all the possible comparisons
     between exons of t1 and introns of t2 and reciprocally (not sym). *)
  let no_exon_over_intron t1 t2 =
    let il1 = intronlist t1 and el1 = t1.exlist and il2 = intronlist t2 and el2 = t2.exlist in
    let lexintrpair1 = ref [] and lexintrpair2 = ref [] in
    let u = List.iter (fun ex1 -> let e1 = (Exon.gbeg ex1,Exon.gend ex1) in List.iter (fun i2 -> (lexintrpair1:=(e1,i2)::(!lexintrpair1))) il2) el1 in
    let u = List.iter (fun ex2 -> let e2 = (Exon.gbeg ex2,Exon.gend ex2) in List.iter (fun i1 -> (lexintrpair2:=(e2,i1)::(!lexintrpair2))) il1) el2 in
    let loverexintrpair1 = List.filter (fun (e,i) -> if (Common.foverlap e i) then true else false) (!lexintrpair1) 
    and loverexintrpair2 = List.filter (fun (e,i) -> if (Common.foverlap e i) then true else false) (!lexintrpair2) in 
      (((List.length loverexintrpair1)==0) && ((List.length loverexintrpair2)==0))
    
  (* note that to be compatible two spliced tr first need to overlap by at least one bp *)
  let compatible_spliced_tr t1 t2 =
    ((overlap t1 t2) && (have_overlapping_and_compatible_introns t1 t2) && (no_exon_over_intron t1 t2))

  (* the three next functions are explicitly used by the comptr algorithm and work on stranded spliced transcript,
     examining their intron list.
     tr_exact says whether t1 is an exact occurence of t2 in terms of intron list
     tr_extension says whether t1 is an extension of t2 in terms of intron list = for all introns of t1 that overlap an intron of t2
     tr_inclusion says whether t1 is included in t2 in terms of intron list
     note that an extension or an inclusion cannot be an exact occurence.
  *)
  let tr_exact t1 t2 =
    let il1 = intronlist t1 and il2 = intronlist t2 in
      (il1 = il2)

  let tr_extension t1 t2 =
      let il1 = intronlist t1 and il2 = intronlist t2 in
      let n1 = List.length il1 and n2 = List.length il2 in
	if (n1 <= n2) then
	  false
	else
	  (il2 = (List.filter (fun i1 -> List.mem i1 il2) il1));;

  let tr_inclusion t1 t2 = 
    let il1 = intronlist t1 and il2 = intronlist t2 in
      ((not (tr_exact t1 t2)) && (List.fold_left (fun bool i1 -> (bool && (List.mem i1 il2))) true il1))


  (* given an output channel, a boolean for ucsc and a transcript, prints the transcript and its exons in gtf format *)
  let print_gtf o u t = 
    if (isnull t) then
      ()
    else
      begin
	(* prints the transcript in gtf format *)
	if u then
	  Printf.fprintf o "%s\t%s\ttranscript\t%i\t%i\t0.0\t%s\t.\tgene_id %s transcript_id %s\n" t.chrom t.cat t.gbeg t.gend (strand_to_string t.str) (value_to_ucsc_value t.gnid) (value_to_ucsc_value t.trid)
	else
	  Printf.fprintf o "%s\t%s\ttranscript\t%i\t%i\t0.0\t%s\t.\tgene_id %s transcript_id %s\n" t.chrom t.cat t.gbeg t.gend (strand_to_string t.str) t.gnid t.trid;
	(* prints the exons in gtf format *)
	List.iter (Exon.print_gtf o u) (t.exlist)
      end
      
  let chrom t = t.chrom
  let gbeg t = t.gbeg
  let gend t = t.gend
  let str t = t.str	
  let cat t = t.cat
  let exlist t = t.exlist
  let trid t = t.trid
  let gnid t = t.gnid

  (* change the exon list of a transcript *)
  let setexlist el t = {t with exlist=el}

  (* changes the coord of the transcript (beg and end) *)
  let setcoord b e t = {t with gbeg=b; gend=e}
end



module Gene = 
struct
  type t = {
	chrom: string;
	gbeg: int;
	gend: int;
	str: strand;
	trarr: Transcript.t array;
	exarr: Exon.t array;
	gnid: string;
  }

  let create chr gb ge st ta ea gn = 
    { chrom=chr; gbeg=gb; gend=ge; str=st; trarr=ta; exarr= ea; gnid=gn}
	  
  let chrom g = g.chrom
  let gbeg g = g.gbeg
  let gend g = g.gend
  let str g = g.str	
  let trarr g = g.trarr
  let exarr g = g.exarr
  let gnid g = g.gnid
end



(* Exon projection. In a gene it is a maximal set of overlapping exons *)
module ExonProj =
struct
  type t = {
	chrom: string;
	gbeg: int;
	gend: int;
	str: strand;
	nbex: int;  (* number of exons if is made of *)
	exarr: Exon.t array;  (* The exons it is made of *)
	noingn: int;    (* number of the exon projection among all the exon projections of the gene (from 5') *)
	begingn: int;   (* begining of the exon projection in the virtual cDNA made by joining all the exon projections of the gene *)
	endingn: int;   (* end of the exon projection in the virtual cDNA made by joining all the exon projections of the gene *)
 }
	  
  let create chr gb ge st ne ea no bgn egn = 
    { chrom=chr; gbeg=gb; gend=ge; str=st; nbex=ne; exarr=ea; noingn=no; begingn=bgn; endingn=egn}
  
  let chrom e = e.chrom
  let gbeg e = e.gbeg
  let gend e = e.gend
  let str e = e.str
  let nbex e = e.nbex
  let exarr e = e.exarr
  let noingn e = e.noingn
  let begingn e = e.begingn 
  let endingn e = e.endingn
end 




(* Segmented (Exon) projection. The different exons of an exon projection define 
   elementary segments called Segmented (Exon) projections.
   The property of a segmented projection is that all its nt are included in the same
   number of transcripts.
*)
module SegProj = 
struct
  type t = {
	chrom: string;
	gbeg: int;
	gend: int;
	str: strand;
	score : float;
	exproj: ExonProj.t;
	coverage: int;   (* number of transcripts it belongs to *)
	noinexproj: int; (* number of the segmented projection in the exon projection from the 5' *)
	begingn: int; (* begining of the segmented projection in the virtual cdna of the gene *)
	endingn: int; (* end of the segmented projection in the virtual cdna of the gene *)
	
 }

  let create chr gb ge st sc ep cov no bgn egn = 
    { chrom=chr; gbeg=gb; gend=ge; str=st; score=sc; exproj=ep; coverage=cov; noinexproj=no; begingn=bgn; endingn=egn}
  
  let setcoverage cov s = {s with coverage=cov}
  let setscore sc s = {s with score=sc}

  let chrom s = s.chrom
  let gbeg s = s.gbeg
  let gend s = s.gend
  let str s = s.str
  let score s = s.score
  let exproj s = s.exproj
  let coverage s = s.coverage
  let noinexproj s = s.noinexproj
  let begingn s = s.begingn
  let endingn s = s.endingn
end


  


(* This function takes as input a segseq of objects that we can associate to a function providing the strand,
   ordered according to strand (forward first, then unknown or unstranded and finally reverse), 
   and simply cuts it wrt strand. Note: this does not modify the input, but returns the result in a new segseq. *)
let cut_according_to_strand_gal strfun s =
  let auxcut a =
    let i = ref 0 and lpos = ref [0] and currstr = ref Forward and lga = Array.length a in 
      while (!i<lga) do
	currstr := (strfun a.(!i));
	while ((!i<lga) && ((strfun a.(!i))=(!currstr))) do
	    incr i;  
	done;
	  if(!i<lga) then
	    lpos := !i::(!lpos);
	done;
      (* Common.print_log (("# I have ")^(string_of_int (List.length (!lpos)))^(" elements in !lpos\n")); *)
      List.rev (!lpos) in 
    SegSeq.setsegment s (auxcut (SegSeq.tank s));;

let cut_according_to_strand_tr s_tr = cut_according_to_strand_gal Transcript.str s_tr;;

let cut_according_to_strand_ex s_ex = cut_according_to_strand_gal Exon.str s_ex;;


(* This function takes as input a SegSeq s containing an array of objects that
   we can associatate to a function providing the chr, and that are ordered
   according to chromosome and cuts the segseq with respect to the chromosome *)
let cut_according_to_chrom_gal chrfun s = 
  let auxcut a =
    let i = ref 0 and lpos = ref [0] and currchr = ref "" and lga = Array.length a in 
      while (!i<lga) do
	currchr := (chrfun a.(!i));
	while ((!i<lga) && ((chrfun a.(!i))=(!currchr))) do
	    incr i;  
	done;
	  if(!i<lga) then
	    lpos := !i::(!lpos);
	done;
      List.rev (!lpos) in 
    SegSeq.setsegment s (auxcut (SegSeq.tank s));;

let cut_according_to_chrom_tr s_tr = cut_according_to_chrom_gal Transcript.chrom s_tr;;

let cut_according_to_chrom_ex s_ex = cut_according_to_chrom_gal Exon.chrom s_ex;;


(* The 2 next functions are auxiliary functions to cut_according_to_position_gal below,
   take as input an array of objects that we can associate to functions providing chromosome,
   beg, end and strand, and cuts them according to position by two different ways: simple 
   (which would be good for objects ordered from 5' to 3'but where the strand is not reverse), 
   and backward (which would be good for objects ordered from 5' to 3' but where the strand IS reverse). 
   Those two functions return the list of positions where to cut in the array of objects a given as input.
   Note that those two functions are inspired from MakeProjection/main.ml (project) which 
   itself was a modified version of DetermineSP/compute_exproj.ml (makeSP). 
   Note: this does not modify the input, but returns the result in a new segseq
*)
let cut_wrtpos_simple_gal chrfun gbegfun gendfun strfun a =
  let i = ref 0 and lpos = ref [0] and currchrom = ref "" and currend = ref 0 and lga = Array.length a in
  while (!i<lga) do
    currchrom := chrfun a.(!i);
    currend := gendfun a.(!i);
    
    while ((!i<lga) && ((chrfun a.(!i))=(!currchrom)) && ((gbegfun a.(!i)) <= !currend)) do
      if ((gendfun a.(!i)) > !currend) then
	currend:= (gendfun a.(!i));
      incr i;
    done;
    if(!i<lga) then
      lpos := !i::(!lpos);
  done;
  List.rev (!lpos)


let cut_wrtpos_backward_gal chrfun gbegfun gendfun strfun a =
  let i = ref 0 and lpos = ref [0] and currchrom = ref "" and currbeg = ref 0 and lga = Array.length a in
    while (!i<lga) do
      currchrom := chrfun a.(!i);
      currbeg := gbegfun a.(!i);
      
      while ((!i<lga) && ((chrfun a.(!i))=(!currchrom)) && ((gendfun a.(!i)) >= !currbeg)) do
	if ((gbegfun a.(!i)) < !currbeg) then
	  currbeg:= (gbegfun a.(!i));
	incr i;
      done;
      if(!i<lga) then
	lpos := !i::(!lpos);
    done;
    List.rev (!lpos)


(* This function takes as input a SegSeq s containing an array of objects that we
   can associated to functions providing chr, beg, end and strand, supposed to be 
   all on the same strand, and according to this strand cuts the segseq with respect 
   to position in a way that depends on the strand (since the objects are ordered 
   from 5' to 3' within a strand here). It modifies the segseq. 
   Note: I should catch the exception for a.(0), not done yet *)
let cut_according_to_position_gal chrfun gbegfun gendfun strfun s =
  let a = SegSeq.tank s in
    try
      (
	match (strfun a.(0)) with
	  | Reverse -> SegSeq.setsegment s (cut_wrtpos_backward_gal chrfun gbegfun gendfun strfun a)
	  | _ -> SegSeq.setsegment s (cut_wrtpos_simple_gal chrfun gbegfun gendfun strfun a)
      )
      with
	| Invalid_argument _ -> SegSeq.make2 [||];; 

let cut_according_to_position_tr s_tr = cut_according_to_position_gal Transcript.chrom Transcript.gbeg Transcript.gend Transcript.str s_tr

let cut_according_to_position_ex s_ex = cut_according_to_position_gal Exon.chrom Exon.gbeg Exon.gend Exon.str s_ex


(* aux function to arrtr_to_mergedtr which makes a string out of list of tr ids 
   that should be merged (comma separated list). Need to catch exception for String.sub 
   in case string is too short *)
let ltrids_to_string ltrids =
  let s= ref "" in
    (* note that here we are adding a space at the end of each line, that is why we do clean_end_string *)
    List.iter (fun trid -> s:=(!s)^(trid)^(",")) ltrids;
    (clean_end_string (String.sub !s 0 ((String.length !s)-1)))


(* takes as input an integer and an array of transcripts that should be merged into one tr
   according to their intron list if spliced and to simple stranded overlap if monoex,
   and outputs a single transcript object with as two attributes the integer as a new id 
   is for this tr and the list of individual transcripts it is made of. 
   Note that the transcripts in the array are supposed to be ordered from 5' to 3' and 
   the exons in each tr are also supposed to be in this order. *)
let arrtr_to_mergedtr i arrtr =
  let ltr = Array.to_list arrtr in
    try
      (
	let t0 = arrtr.(0) in
	let lallex_ordered = List.fold_left (fun lcumul lex -> List.merge Exon.compare_pos lcumul lex) [] (List.map Transcript.exlist ltr) in
	let segseq_allex = SegSeq.make2 (Array.of_list lallex_ordered) in
	let segseq_allex_cut = cut_according_to_position_ex segseq_allex in 
	let larrex = SegSeq.elements segseq_allex_cut in
	let lprojex = List.map exarr_to_projex larrex in
	(* transcript id will be the new tr id and ene id can be a place for remembering which indiv tr are merged here  *)
	let trid=string_of_int i and gnid=ltrids_to_string (List.map Transcript.trid ltr) in
	  Transcript.create 
	    (Transcript.chrom t0)
	    (List.fold_left min (Transcript.gbeg t0) (List.map Transcript.gbeg ltr))
	    (List.fold_left max (Transcript.gend t0) (List.map Transcript.gend ltr))
	    (Transcript.str t0)
	    (Transcript.cat t0)
	    (List.map (Exon.settrgn trid gnid) lprojex)             
	    trid 
	    gnid
      )
    with
      | Invalid_argument _ -> Transcript.null (* means the array of tr is empty, in this case we create a null object *)
  

let listtr_to_mergedtr i ltr =
  try
    (
      let t0 = List.hd ltr in
      let lallex_ordered = List.fold_left (fun lcumul lex -> List.merge Exon.compare_pos lcumul lex) [] (List.map Transcript.exlist ltr) in
      let segseq_allex = SegSeq.make2 (Array.of_list lallex_ordered) in
      let segseq_allex_cut = cut_according_to_position_ex segseq_allex in 
      let larrex = SegSeq.elements segseq_allex_cut in
      let lprojex = List.map exarr_to_projex larrex in
	(* transcript id will be the new tr id and ene id can be a place for remembering which indiv tr are merged here  *)
      let trid=string_of_int i and gnid=ltrids_to_string (List.map Transcript.trid ltr) in
	Transcript.create 
	  (Transcript.chrom t0)
	  (List.fold_left min (Transcript.gbeg t0) (List.map Transcript.gbeg ltr))
	  (List.fold_left max (Transcript.gend t0) (List.map Transcript.gend ltr))
	  (Transcript.str t0)
	  (Transcript.cat t0)
	  (List.map (Exon.settrgn trid gnid) lprojex)             
	  trid 
	  gnid
    )
  with
    | Failure _ -> Transcript.null (* means the list of tr is empty, in this case we create a null object *)
  
