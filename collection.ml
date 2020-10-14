(* collection.ml : define the SegSeq object (in the form of a module) = a segmented sequence, i.e. a set of objects
   of the same type that are ordered and that are segmented into several ones.
   For example several genes along the sequence of a given genome *)


module SegSeq = 
struct
  type 'a t = { 
    tank : 'a array ; (* array to store the set of ordered objects *)  
    ord : int list    (* segmentation positions within the previous array *)
  }

  let tank s = s.tank
  let ord s = s.ord
  let make size v = { tank= Array.make size v ; ord = []}
  let make2 t = {tank=t; ord = []}
  let make3 t o = {tank=t; ord = o}
  let create = make 

  let isnull s = (s.tank = [||])

  let hd s = 
    match s.ord with 
      | [] -> failwith "SegSeq.hd"
      | i::_ -> i
	  
  exception Seg of int 
    

  let setsegment s l = {s with ord=l}

  (* further partition a segseq according to a predicate p *)
  let partition p s = 
    let rec apart deb fin = 
      try 
	for i = deb to fin do 
	  if not (p s.tank i) then raise (Seg(i))  (* p takes an array and the current position as input *)
	done;
	[]
      with Seg i -> i::(apart (i+1) fin)
 
    in
    let rec lpart = function
      |[] -> [] 
      |[i]-> [apart i ((Array.length s.tank) -1)]
      |i::(j::r as l) -> (apart i j)::(lpart l)
    in
      setsegment s (List.flatten (lpart s.ord))

 
  (* returns the intervals corresponding to the segments of the segseq *)
  let intervalles s =
    let rec aux l n =  
      match l with
	|[] -> [];
	|[t] -> [(t,n)];
	|t1::t2::q -> (t1,t2)::(aux (t2::q) n) in
      aux s.ord ((Array.length s.tank)-1)


  let single a = {tank = a ; ord =[]}

 (*  let fold f s e =  *) 


  (* Converts a segseq into a list of arrays. 
     Made tail recursive in order not to crash for big tanks.
  *)
  let elements s = 
    let rec lambda acc l = match l with
      |[] -> []
      |[i]-> List.rev ((Array.sub s.tank i ((Array.length s.tank)-i))::acc)
      |i::(j::_ as q)-> lambda ((Array.sub s.tank i (j-i))::acc) q  
    in
      lambda [] (s.ord)
 
  (*  Old version, but be careful makes a stack overflow (in the form of a segmentation fault) for very long segseq
      export OCAMLRUNPARAM='l=100000000000000000000000000000000000000000k'
      does not change anything (wrong syntax here?)
      let elements s = 
      let rec lambda = function 
      |[] -> []
      |[i]-> [Array.sub s.tank i ((Array.length s.tank)-i)]
      |i::(j::_ as q)-> (Array.sub s.tank i (j-i))::(lambda q)  
      in
      lambda (s.ord) 
  *)  

  (* f is a function that takes an element of the type of what is in tank *)
  let map f s =
    make3 (Array.map f s.tank) s.ord


  (* converts a list of arrays into a segseq *) 
  let tosegseq lt =
    let i = ref 0 and s = ref 0 and n = List.length lt in
    let tord = Array.create n 0 in
      
      while(!i<n) do
	tord.(!i) <- !s;
	s := !s + Array.length (List.nth lt !i);
	incr i;
      done;
      
      setsegment (make2 (Array.of_list (List.flatten (List.map Array.to_list lt)))) (Array.to_list tord)
	
end 
    

    

      
