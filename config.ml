(* config.ml *)
(* module to manage the parameters input to the program, and to display the help and usage *)

   

(***********************************************************)
(* Structure for the parameters of the trmerge program     *)
(***********************************************************)
type 'a context_t =
	{ 
	  mutable file1: string;	        (* input file1 = exon gtf file of transcripts to assess *)
 	  mutable file2: string;	        (* input file2 = exon gtf file of reference transcripts *)
 	  mutable maxnbseq: int;                (* maximum number of genomic sequences in the input files *)
	  mutable outfile: string;              (* output file = a three column tsv file with transcript id of file1, its class wrt file2, and a list of tr of file2 *)
	  mutable verbose: bool;                (* if we want the output to be displayed on stdout *)
	  mutable sorted: bool;                 (* boolean for whether file is already sorted according to strand, chr, start, end 
						   (sort -k7,7 -k1,1 -k4,4n -k5,5n) *)
	} 



(********************************************************)
(* trmerge context, these are the default parameters    *)
(********************************************************)
let context = 
  {	
    file1 = "";
    file2 = "";
    maxnbseq=50;
    outfile="";
    verbose=false;
    sorted=false;
  };;


let usage = 
" 
                     ***********************************************          
                     *   comptr - version v1.1 (February 2014)     *
                     *               Sarah Djebali                 *
                     ***********************************************

Usage: "^(Filename.basename (Sys.argv.(0)))^" file1 file2 [options] 

Takes as input two gff version 2 files, the first one with transcript structures to assess and the second 
one with reference transcript structures against which to do this assessment (only exon rows will be considered
but they have to contain transcript_id in column 11 and its value in column 12 otherwise it will exit). 
This program will only assess stranded spliced transcripts coming from the first file. It will output 
a three column tsv file where the first column is the assessed transcript id, the second column is its 
class with respect to file2 transcripts, and the third column is a comma separated list of file2 transcript ids 
corresponding to the assessed transcript.
The output classes are the following and are chosen in the following order in case of multiple possible assignments:
- Unstranded: the assessed transcript is unstranded,
- Monoexonic: the assessed transcript is mono-exonic,
- Intergenic_or_antisense: the assessed transcript is stranded and spliced but does not strandedly overlap any reference transcript,
- Exact: there is a reference transcript with all its introns equal to the assessed transcript and reciprocally,
- Inclusion: there is a reference transcript with a subset of introns corresponding to all introns of the assessed transcript,
- Extension: there is a reference transcript with all its introns equal to the assessed transcript but the assessed
  transcript has additional introns,
- Overlap: there is a reference transcript overlapped by the assessed transcript,

** file1 and file2 must be provided in gff version 2 format and should contain exon rows (with transcript_id in field no 12).
** [options] can be:
   -o outfile:    outfile is the name of the gff file the user wants the output to be provided in. 
                  -> default is file1_comp.tsv.

   -v:            provides the output result in the standard output rather than in an output file.
                  -> default is unset.
   
   -s nbseq:      nbseq is an upper bound for the number of sequences you have in your input gff files.
                  -> default is 50.

   -so:           comptr does not require the input file to be sorted, however this option enables to skip
                  the input file sorting in case this file is already sorted according to strand, chromosome, start, end. 
                  Note: this sorting could be performed outside compmerge using the unix sort command: 
                  sort -k7,7 -k1,1 -k4,4n -k5,5n file > sortedfile
** Please report any bug to sarahqd@gmail.com        
"




(***********************************************************************)
(* Read the arguments from the command line and updates the context    *)
(***********************************************************************)
let read_commandline () =
  let u = try 
      begin
	context.file1 <- Sys.argv.(1);
	context.file2 <- Sys.argv.(2);
      end
    with
      | Invalid_argument s -> Common.print_error usage
  in
  
  (* we start reading the arguments from the 3rd one since the two first are compulsory and are the input files *)
  let argnum = ref 2 and ok = ref true in

  (* This function returns the next argument. mustbe says whether this argument has to exist.
     In general mustbe is true since the arguments go by pairs *)
  let getarg mustbe =
    incr argnum; 
    if !argnum < (Array.length Sys.argv) then 
      Sys.argv.(!argnum)
    else if mustbe then 
      raise Not_found 
    else
      ""
  in
    (* Actually reading each of the arguments *)
    try 
      while (!ok) do
	match (getarg false) with
	  | "-o"  -> context.outfile <- getarg true
	  | "-v"  -> context.verbose <- true
	  | "-s"  -> context.maxnbseq <- int_of_string (getarg true)
	  | "-so" -> context.sorted <- true
	  | "-h"
	  | "--help" -> Printf.printf "%s\n" usage; exit 1; 
	  | ""	-> ok := false
	  | s	-> failwith s
      done;
      Common.print_log "# Command line has been read\n";
    with
      | Not_found -> Common.print_error ("Missing parameter at the end of the command line\n"^usage^"\n");
      | Failure "int_of_string" -> Common.print_error ("Syntax error (incorrect integer value)\n"^usage^"\n");
      | Failure s -> Common.print_error ("Syntax error ("^s^")\n"^usage^"\n");;


