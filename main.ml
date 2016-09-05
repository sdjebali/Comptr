(* main.ml *)

open Common
open Collection
open Config
open Feature
open Input


(* the class of an assessed transcript compared to a set of reference transcripts
   this is the result of the comparison algorithm comptr is trying to achieve
*)
type trclass = | Unstranded
    | Monoexonic
    | IntergenicAntisense
    | Exact
    | Inclusion
    | Extension
    | Overlap;;


let class_to_string cl = 
  match cl with
    | Unstranded -> "Unstranded"
    | Monoexonic -> "Monoexonic"
    | IntergenicAntisense -> "Intergenic_or_antisense"
    | Exact -> "Exact"
    | Inclusion -> "Inclusion"
    | Extension -> "Extension"
    | Overlap -> "Overlap";;
    

(* 
   The function overlap takes as input two list of transcripts from the same strand and chromosome,
   and ordered accprding to coordinates (start, end), and outputs a list of pairs of the same size as the first 
   but with pairs where the first element is the transcript in list1 and the second element is the list 
   of tr of list2 that overlap it. It goes over the two lists at the same time to be faster.  
   if trlist1 is empty then the list returned is also empty
   if trlist2 is empty then the list returned has the same length as trlist1 but with empty list in second element of the pairs
   need to check that this algo is correct on a simple example in the interpreter
   *)                                                                                                      
let overlap trlist1 trlist2 =                                         
  let nbtr1 = List.length trlist1 and i2 = ref 0 and nbtr2 = List.length trlist2 and pairlist1 = ref [] in
    for i1=0 to (nbtr1-1) do                                                                                 
      let tr2overtr1 = ref [] in                                                                      
        while (((!i2)<nbtr2) && ((Transcript.gend (List.nth trlist2 (!i2)))<(Transcript.gbeg (List.nth trlist1 i1)))) do        
          incr i2;
        done;
        let firsti2 = !i2 in
          if (((!i2)<nbtr2)&&((Transcript.gbeg (List.nth trlist2 (!i2)))<=(Transcript.gend (List.nth trlist1 i1)))) then
            begin
              while ((!i2)>=0&&(Common.foverlap (Transcript.interval (List.nth trlist1 i1)) (Transcript.interval (List.nth trlist2 (!i2))))) do
                decr i2;
              done;
              incr i2;
              while (((!i2)<nbtr2)&&((Transcript.gbeg (List.nth trlist2 (!i2)))<=(Transcript.gend (List.nth trlist1 i1)))) do
                if (Common.foverlap (Transcript.interval (List.nth trlist1 i1)) (Transcript.interval (List.nth trlist2 (!i2)))) then
                    tr2overtr1:=(((List.nth trlist2 (!i2)))::(!tr2overtr1));
                incr i2;
              done; (* end of while (((!i2)<nbtr2)&&(Transcript.gbeg (trlist2.(!i2)))<=(Transcript.gend (trlist1.(i1)))) *)
              i2:=firsti2;
            end;  (* end of if (((!i2)<nbtr2)&&(Transcript.gbeg (trlist2.(!i2)))<=(Transcript.gend (trlist1.(i1)))) *)
          pairlist1:=((List.nth trlist1 i1),(List.rev (!tr2overtr1)))::(!pairlist1);
    done;  (* end of for i1=0 to (nbtr1-1) *)
    List.rev (!pairlist1);;


(* by definition the two lists are not empty since come from arrays
   given by SegSeq.elements *)
let find_trlist_with_same_strand_and_chrom trll trlist = 
  List.find (fun trl -> Transcript.same_strand_and_chrom (List.hd trlist) (List.hd trl)) trll;;
      

(* Given a transcript tr1 and a list of transcripts trlist2 that overlap it on the same strand 
   provide a class for tr1 as well as the sublist of tr from trlist2 associated to it.
   Priority is give to unstranded tr, then monoex, then in case of overlap priority is given
   to exact, then inclusion and then extension. I will associate to tr1 its class and its list 
   of tr using a triplet as output of this function.
*)
let decide_tr_class (tr1,listtr2) =
  let ltrexact = ref [] and ltrinclusion = ref [] and ltrextension = ref [] in
  if (((Transcript.str tr1)!=Forward) && ((Transcript.str tr1)!=Reverse)) then
    (tr1,Unstranded,[])
  else
    begin
      if ((List.length (Transcript.exlist tr1))==1) then
	(tr1,Monoexonic,[])
      else
	begin
	  if ((List.length listtr2)==0) then
	    (tr1,IntergenicAntisense,[])
	  else
	    begin
	      ltrexact:=(List.filter (Transcript.tr_exact tr1) listtr2);
	      if ((List.length (!ltrexact)) != 0) then
		(tr1,Exact,(!ltrexact))
	      else
		begin
		  ltrinclusion:=(List.filter (Transcript.tr_inclusion tr1) listtr2);
		  if ((List.length (!ltrinclusion)) != 0) then
		    (tr1,Inclusion,(!ltrinclusion))
		  else
		    begin
		      ltrextension:=(List.filter (Transcript.tr_extension tr1) listtr2);
		      if ((List.length (!ltrextension)) != 0) then
		        (tr1,Extension,(!ltrextension))
		      else
			(tr1,Overlap,listtr2)
		    end
		end
	    end
	end
    end;;


let print_triplet_tr1_class_tr2list_as_tsv_row o (tr1,cl,listtr2) =
  Printf.fprintf o "%s\t%s\t%s\n" (Transcript.trid tr1) (class_to_string cl) (list_to_string2 (List.map Transcript.trid listtr2));;


(**********************************************  MAIN FUNCTION *************************************)

(* this is the function that executes the outest actions *)
let comptr_main () =
  (* Read the arguments on the command line *)
  read_commandline ();
  let u = Common.print_log (("# Input file1 is ")^(context.file1)^("\n")) in
  let u = Common.print_log (("# Input file2 is ")^(context.file2)^("\n")) in
  
  (* Sort input files according to strand and then gx coordinates (chr, start, end), and put the result in  
     an intermediate file tmp_sorted_file placed by default in the /tmp directory since a /tmp directory always 
     exists on all systems, but it may be interesting to add an option to enable the user
     to specify his own tmp directory (in case of huge file for example). This argument
     would then need to be passed to the make_temp_file_name function *)
  let tmp_sorted_file1 = Common.make_sorted_temp_file_if_user_wants_stranded true context.sorted context.file1 in
  let tmp_sorted_file2 = Common.make_sorted_temp_file_if_user_wants_stranded true context.sorted context.file2 in
  let u= Common.print_log ("# I have treated (sorted and put in temp file) the input files according to what the user wants\n") in

  (* Open input channel *)
  let inchan1 = open_in tmp_sorted_file1 and inchan2 = open_in tmp_sorted_file2 in
  
  (* Open output channel *)
  let outchan = 
    if (context.verbose) then
      stdout
    else
      begin
	try
	  (open_out context.outfile) 
	with
	  | Sys_error s -> open_out ((context.file1)^("_comp.tsv"))
      end
  in
 
  (* 1. Read from input channel and make list of Transcript objects containing themselves exon objects in their exon lists.
     note: those tr are already ordered according to strand first and then from 5' to 3' and have their exons already 
     ordered from 5' to 3'. *)
  let trlist1 = read_gff_into_tr_list2 inchan1 and trlist2 = read_gff_into_tr_list2 inchan2 in
  let u = Common.print_log (("# I have ")^(string_of_int (List.length trlist1))^(" initial transcripts in file1\n")) in
  let u = Common.print_log (("# I have ")^(string_of_int (List.length trlist2))^(" initial transcripts in file2\n")) in
  
  (* 2. Divide each list into sublists using segseq and cut according to strand first and then according to chromosome
     since the overlap algorithm works for a given strand and chromosome *) 
  let segseq1 = SegSeq.make2 (Array.of_list trlist1) and segseq2 = SegSeq.make2 (Array.of_list trlist2) in
  let segseq1_cut_by_strand = cut_according_to_strand_tr segseq1 and segseq2_cut_by_strand = cut_according_to_strand_tr segseq2 in
  let larr1_cut_by_strand = SegSeq.elements segseq1_cut_by_strand and larr2_cut_by_strand = SegSeq.elements segseq2_cut_by_strand in  (* list of arrays where each array
																	 corresponds to a strand *)
  let u = Common.print_log ("# I have divided file1 and file2 transcripts intro sublists according to strand\n") in

  let lsegseq1_cut_by_strand_and_chr = List.map (fun a -> cut_according_to_chrom_tr (SegSeq.make2 a)) larr1_cut_by_strand in 
  let lsegseq2_cut_by_strand_and_chr = List.map (fun a -> cut_according_to_chrom_tr (SegSeq.make2 a)) larr2_cut_by_strand in  (* for each strand we have a segseq that
																 we further cut according to chr *)
  let llarr1_cut_by_strand_and_chr = List.map SegSeq.elements lsegseq1_cut_by_strand_and_chr in
  let llarr2_cut_by_strand_and_chr = List.map SegSeq.elements lsegseq2_cut_by_strand_and_chr in (* list of list of arrays for cut according to strand and then chr, For each strand,
												   we have a list of arrays where each array is an array of tr of the same chr *)  

  let larr1_cut_by_strand_and_chr = List.flatten llarr1_cut_by_strand_and_chr in
  let larr2_cut_by_strand_and_chr = List.flatten llarr2_cut_by_strand_and_chr in    (* only one list of arrays, where each array correspond to a different combination of 
										       strand and chromosome *)
   
  let llist1_cut_by_strand_and_chr = List.map Array.to_list larr1_cut_by_strand_and_chr in
  let llist2_cut_by_strand_and_chr = List.map Array.to_list larr2_cut_by_strand_and_chr in   (* list of list instead of list of arrays *)
  let u = Common.print_log ("# I have divided file1 and file2 transcripts intro sublists according to strand and chromosome\n") in


  (* 3. For each strand and chromosome, associate to each transcript of the first set the list of transcripts
     from the second set of the same strand and chromosome that overlap it. *)
  let llist1_of_pairs = List.map (fun trlist1 -> overlap trlist1 (try (find_trlist_with_same_strand_and_chrom llist2_cut_by_strand_and_chr trlist1) with Not_found -> [])) llist1_cut_by_strand_and_chr in
  let list1_of_pairs = List.flatten llist1_of_pairs in
  let u = Common.print_log ("# For each strand and chromosome I have associated each transcript of file1 to a set of strandedly overlapping transcripts from file2\n") in


  (* 4. By comparing the transcript of list1 to all the tr of list2 that overlap it, decide on its class *)
  let list1_of_tr_with_class_and_tr2 = List.map decide_tr_class list1_of_pairs in
  let u = Common.print_log ("# For each transcript of file1, I have decided on its class and assiated transcripts from file2\n") in

  
  (* 5. Print the tr of set 1 with its class and associated tr of set1 as a 3 column tsv file *)
  let u = List.iter (print_triplet_tr1_class_tr2list_as_tsv_row outchan) list1_of_tr_with_class_and_tr2 in
  let u = Common.print_log ("# I have printed the transcripts of file1 with their class and associated transcripts from file2, in a 3 column tsv file\n") in

  let u = Common.remove_sorted_temp_file_if_user_wants context.sorted tmp_sorted_file1 in
  let u = Common.remove_sorted_temp_file_if_user_wants context.sorted tmp_sorted_file2 in
  let u = Common.print_log ("# I have removed the temporary sorted files in case they have been created\n") in
  let u = Common.print_log ("# comptr did its work !\n") in
    flush outchan;;   (* note: tr of file1 are printed in the reverse order but is it a problem? *)


comptr_main ();;

 
