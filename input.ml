(* input.ml *)
(* reads input files and store the data in appropriate structures *)


open Common
open Feature



(* record_of_line_gff_flex takes as input a line from a gff file as a list of strings 
   split by the "\t" separator (flexible gff), and outputs a feature object.*)
let record_of_line_gff_flex lline =
  let tab = Array.of_list lline in
  let n = Array.length tab in
    if (n <= 7) then
      begin
	Common.print_error "# record_of_line_gff_flex: syntax error in your gff file: it should contain at least 8 fields separated by tabs"; 
	Feature.null
      end
    else
      begin
	if (n==8) then
	  begin
	    Feature.create 
	      tab.(0) 
	      (int_of_string tab.(3)) 
	      (int_of_string tab.(4)) 
	      (strand_of_string tab.(6)) 
	      tab.(2) 
	      tab.(1) 
	      tab.(5)
	      tab.(7)
	      []
	  end
	else
	  begin
	    let attl = split ' ' (Common.clean_beg_string (Common.clean_end_string tab.(8))) in
	    let natt = List.length attl in
	      if((natt mod 2)!=0) then
		begin
		  Common.print_error "# record_of_line_gff_flex: syntax error in your gff file: it should contain at least 8 fields separated by tabs, followed by a list of attibutes of the form [key space value] separated by spaces. There should not be any space at the end"; 
		  Feature.null
		end
	      else
		begin
		  let rec toattlist = function
		      [] -> []
		    | k::v::q -> (k,v)::(toattlist q)
		  in
		    (* try( *)
		      Feature.create 
		      tab.(0) 
		      (int_of_string tab.(3)) 
		      (int_of_string tab.(4)) 
		      (strand_of_string tab.(6)) 
		      tab.(2) 
		      tab.(1) 
		      tab.(5)
		      tab.(7)
		      (toattlist attl)
		    (* )with
		    //  | Failure _ -> let s=((tab.(0))^("\t")^(tab.(1))^("\t")^(tab.(2))^("\t")^(tab.(3))^("\t")^(tab.(4))^("\t")^(tab.(5))^("\t")^(tab.(6))^("\t")^(tab.(7))^("\t")^(tab.(8))) in Common.print_pb (("# ")^s^("\n")); Feature.null *)
		end
	  end
      end


(* This functions reads from an input channel and makes a list of transcript objects
   from the exon lines seen in the channel. Note: it requires that exon lines contain 
   transcript id in field no 12, preceeded by transcript_id in field no 11 otherwise 
   it will die. The exons read will be stored as exon objects in
   the list of exons a tr object contains (get inspired from make_exon_list in input.ml 
   of makesp code). The maximum number of transcripts is 1 million here but should
   be taken outside. The hashtable has tr ids for keys and tr as values. 
   Also it returns a list of tr which have correct beg, end and exon list, where
   the exon list of each tr is ordered from 5' to 3', and a list of tr which are themselves
   ordered from 5' to 3' *)
let read_gff_into_tr_list inchan =
  let ltr= ref [] and htr = Hashtbl.create 1000000 and stop = ref false in
    while (not !stop) do
      (
	try
	  let currline = (split '\t' (input_line inchan)) in
	  (* let u = Common.print_log (("# Current line has ")^(string_of_int (List.length currline))^(" fields\n")) in *)
	    if ((List.nth currline 2)="exon") then
	      begin
		let currfeat= record_of_line_gff_flex currline in
		let currex = feature_to_exon currfeat in
		  if (Hashtbl.mem htr (Exon.trid currex)) then
		    begin
		      let t = Hashtbl.find htr (Exon.trid currex) in
		      let lex = Common.insert_right_place false currex (Transcript.exlist t) Exon.compare_pos in
		      let gbeg = List.fold_left min 3000000000 (List.map Exon.gbeg lex) and gend = List.fold_left max (-1) (List.map Exon.gend lex) in
			Hashtbl.replace htr (Exon.trid currex) (Transcript.setcoord gbeg gend (Transcript.setexlist lex t)) 
		    end
		  else
		    Hashtbl.add htr (Exon.trid currex) (Transcript.create (Exon.chrom currex) (Exon.gbeg currex) (Exon.gend currex) (Exon.str currex) (Exon.cat currex) [currex] (Exon.trid currex) (Exon.gnid currex))
	      end
	with
	  | Not_found -> Common.print_error "# At least one exon row does not have transcript_id in field no 11!"
	  | Failure _ -> Common.print_error "# Incorrectly formatted gtf file with not even 3 fields separated by tabs!"
	  | End_of_file -> stop:=true
      )
    done;
    (* output tr list here, from the hashtable, with the tr 1) having their exons ordered from 5' to 3', 2) with correct start and end 
       3) all tr in the ltr list being ordered according to strand first and then from 5' to 3' *) 
    let u = Common.print_log (("# Hashtbl of transcripts has ")^(string_of_int (Hashtbl.length htr))^(" keys/transcripts\n")) in
      Hashtbl.iter (fun k v -> (ltr:=(v::(!ltr)))) htr;   
      (* Transcript.print_gtf stdout true v; ltr:=(Common.insert_right_place false v (!ltr) (Transcript.compare_strand_and_pos)) *)
      (* List.iter (Transcript.print_gtf stdout true) (!ltr);  *)
      List.stable_sort (Transcript.compare_strand_and_pos) (!ltr);;    
    
 
(* This is exactly the same function as read_gff_into_tr_list except that the transcripts are ordered
   by genomic order and not from 5' to 3' in the resulting list. Also their exons are ordered by
   genomic order and not from 5' to 3' as was used in the compmerge program. *)
let read_gff_into_tr_list2 inchan =
  let ltr= ref [] and htr = Hashtbl.create 1000000 and stop = ref false in
    while (not !stop) do
      (
	try
	  let currline = (split '\t' (input_line inchan)) in
	  (* let u = Common.print_log (("# Current line has ")^(string_of_int (List.length currline))^(" fields\n")) in *)
	    if ((List.nth currline 2)="exon") then
	      begin
		let currfeat= record_of_line_gff_flex currline in
		let currex = feature_to_exon currfeat in
		  if (Hashtbl.mem htr (Exon.trid currex)) then
		    begin
		      let t = Hashtbl.find htr (Exon.trid currex) in
		      let lex = Common.insert_right_place false currex (Transcript.exlist t) Exon.compare_pos_gx in
		      let gbeg = List.fold_left min 3000000000 (List.map Exon.gbeg lex) and gend = List.fold_left max (-1) (List.map Exon.gend lex) in
			Hashtbl.replace htr (Exon.trid currex) (Transcript.setcoord gbeg gend (Transcript.setexlist lex t)) 
		    end
		  else
		    Hashtbl.add htr (Exon.trid currex) (Transcript.create (Exon.chrom currex) (Exon.gbeg currex) (Exon.gend currex) (Exon.str currex) (Exon.cat currex) [currex] (Exon.trid currex) (Exon.gnid currex))
	      end
	with
	  | Not_found -> Common.print_error "# At least one exon row does not have transcript_id in field no 11!"
	  | Failure _ -> Common.print_error "# Incorrectly formatted gtf file with not even 3 fields separated by tabs!"
	  | End_of_file -> stop:=true
      )
    done;
    (* output tr list here, from the hashtable, with the tr 1) having their exons in genomic order, 2) with correct gbeg and gend and
       3) all tr in the ltr list being ordered according to strand first and then genomically *) 
    let u = Common.print_log (("# Hashtbl of transcripts has ")^(string_of_int (Hashtbl.length htr))^(" keys/transcripts\n")) in
      Hashtbl.iter (fun k v -> (ltr:=(v::(!ltr)))) htr;   
      (* Transcript.print_gtf stdout true v; ltr:=(Common.insert_right_place false v (!ltr) (Transcript.compare_strand_and_pos)) *)
      (* List.iter (Transcript.print_gtf stdout true) (!ltr);  *)
      List.stable_sort (Transcript.compare_strand_and_pos_gx) (!ltr);;    
    
 

(* The function make_seq_nbfeat_hash reads a gffflex file and outputs a hashtable
   which keys are the sequence name of the features read, and which values are the
   number of features associated to this sequence in the file *)
let make_seq_nbfeat_hash maxseq infile =
  (* create the hashtable *)
    let h = Hashtbl.create maxseq and inchan = open_in infile and stop = ref false in
      while (not !stop) do
	(
	  try
	    let seq = List.hd (split '\t' (input_line inchan)) in
	      if (not (Hashtbl.mem h seq)) then
		Hashtbl.add h seq (ref 1)
	      else
		incr (Hashtbl.find h seq)
	  with
	    | End_of_file -> stop:=true
	)
      done;
      h
      


(* The function fill_feat_arr_hash takes as input a gffflex file and a hashtable
   corresponding to the sequences and the associated number of features in the gffflex file
   and fills in the hastable with the correct features *)
let fill_feat_arr_hash h infile =
  let inchan = open_in infile and stop = ref false and hintref = Hashtbl.create (Hashtbl.length h) and a= Array.make 10 Feature.null in
    while (not !stop) do
	(
	  try
	    let il = input_line inchan in
	    let currline = (split '\t' il) in
	    let seq = List.hd currline in
	      begin
		if (not (Hashtbl.mem hintref seq)) then
		  Hashtbl.add hintref seq (ref 0)
		else
		  incr (Hashtbl.find hintref seq);
		(try 
		    (Hashtbl.find h seq).(!(Hashtbl.find hintref seq)) <- (record_of_line_gff_flex currline)
		  with 
		    | Invalid_argument s -> Common.print_error (("# Index out of bound problem in (Hashtbl.find h seq).(!(Hashtbl.find hintref seq)) of fill_feat_arr_hash"))
		)
	      end
	  with
	    | End_of_file -> stop:=true
	    | Failure "int_of_string" -> Common.print_error (("# Conversion problem from string to integer in fill_feat_arr_hash"))
	)
      done;
      h


(* This function reads the lines of a gff file and outputs a hashtbl of arrays of features.
   This function takes the name of the input file and not the input channel, 
   since it needs to first remember the name of each sequence and its associated
   number of input features. *)
let make_feat_arr_hash maxseq infile = 
  (* Retrieve the names of the sequences and the number of features associated to them reading the input file 
     and put this info in a first preliminary hashtable *)
  let hseq = make_seq_nbfeat_hash maxseq infile in
  let u = Common.print_log "# Names of sequences retrieved\n" in
    
  (* Create the actual hashtable of feature arrays *)
  let h = Hashtbl.create (Hashtbl.length hseq) in
  let u = Common.print_log "# Hashtable of feature arrays created\n" in
    
    (* Fill the actual hashtbl with feature arrays of correct size *)
    Hashtbl.iter (fun k nbref -> Hashtbl.add h k (Array.make (!nbref) Feature.null)) hseq;
    Common.print_log "# Hashtable of feature arrays filling - 1st step\n";

    (* Fill the arrays with the correct features re-reading the input file *)
    fill_feat_arr_hash h infile;
    Common.print_log "# Hashtable of feature arrays filled\n";

    (* Sort the actual feature arrays *)
    Hashtbl.iter (fun k a -> Array.sort Feature.compare_coord a) h;
    Common.print_log "# Hashtable of feature arrays sorted\n";
    h;;
    

