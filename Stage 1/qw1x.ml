module R = Random
module P = Printf

let allow_prints = false ;;
let loi_n1 x coeff = let n = x +. coeff in (*let n = x /. coeff in 	distance adimensionnée*)					
					 (n /. (n +. 1.));;


let proba : float -> bool = fun p -> R.float 1. < p;;(* renvoie true avec probabilité p*)
type ptcl =	{x : int; y : int; dx : int; dy : int}


let timeflow : int -> int -> ptcl ref -> int list ref -> int list ref -> int list ref -> unit = fun t xd q sptm coll txplot ->
	let p = !q in let xx = p.x + p.dx in 
	q := {x=xx; y= p.y + p.dy; dx=p.dx; dy=p.dy};
	let collision = ref false in
	
	sptm := 0 :: (List.filter (fun x -> x <> 0 && (*s'il y a eu une collision*)
	(let un_metre = float_of_int @@ xd  in (*un mètre -> coefficient à changer pour modifier l'intensité de la gravité*)
	let n = float_of_int x in (*distance dimensionnée, mesurée en epsilons*)
	proba @@  un_metre *. (1. -. (loi_n1 n un_metre))  ))

    @@ List.map (fun x -> if xx = x then (collision := true; 0) else (x + 1)) !sptm);
	
	if !collision then(*on a eu collision avec un ralentisseur*)
		begin coll := xx :: !coll; (*Printf.fprintf stdout "Collision n°%d au temps %d en position %d! \n" (List.length !coll) t p.x;*) q := {x=xx - 2* p.dx; y=p.y+p.dy; dx=p.dx ; dy=p.dy}; end
	else if (!q).x < 2 then (txplot := 1 :: !txplot ; raise Exit);
	
	txplot := (!q).x :: !txplot 		
;;

let init_p x0 y0 dx0 dy0 = {x=x0;y=y0;dx=dx0;dy=dy0}

let print_list : int list -> unit = fun l -> print_newline (); List.iter (fun x -> print_int x; print_string "|") l ;;

let rec entassement x dv = function
	|[] -> [dv,x]
	|s :: q when s = x -> entassement x (dv + 1) q
	|t :: q -> (dv,x) :: (entassement t 1 q)
;;


let pyplot : int list -> unit = fun l -> ( print_newline (); List.iter (fun c -> let x,y = c in print_int y; for i = 1 to x/100 do
print_string "█" done; print_newline ()) (entassement 0 0 (List.map (fun t-> t/10) l))
)
;;


let main : int -> int -> (int list * (int list)) = fun steps x_dep -> 
R.self_init ();
	if allow_prints then print_float (loi_n1 1. 1.);
	let sptm = ref [] in
	let p = ref @@ init_p x_dep 0 (-1) 0 in

	let coll = ref [] in 
	let txplot = ref [] in
 	try
	for t = 0 to steps - 1 do
		timeflow t x_dep p sptm coll txplot
	done
	;
	
	if allow_prints then
	(print_newline ();
	pyplot !coll);
	
	print_newline ();
	(*print_list !coll;*)
	Printf.fprintf stderr "Le trou Noir nous a repoussé avec %d collisions. On finit en position %d\n" (List.length !coll) (!p).x;
	(!coll,!txplot)
	with Exit -> 
		if allow_prints then (print_newline (); pyplot !coll);
		Printf.fprintf stderr "Le trou Noir a été atteint en %d collisions\n" (List.length !coll); (!coll,!txplot) 
 ;;

(* module graphique de python à utiliser*)
let plot_coll : int -> out_channel -> (int list * (int list)) -> unit = fun k outc (coll, xplot) -> 
	let ak = "a" ^ string_of_int k in
	let str1 = ak ^ " = np.array([], dtype = int)" in
	let str2 = (ak ^ " = np.array([") in
	let str3 = ("])\n#graph(" ^ ak ^ ")\n") in 
	match List.rev xplot with
		|[] -> P.fprintf outc "%s" str1
		|t :: q -> 
			(P.fprintf outc "%s"  str2; P.fprintf outc "%d" t; 
			List.iter (fun x -> P.fprintf outc "," ;P.fprintf outc "%d" x) q; P.fprintf outc "%s" str3 )(*\ngeodesics(a[0]+1,len(a))*)
in


let m = 225 in(*5000 is good - 5000000 is fine though*)
let n = 500 in (*  2 times m is good *)
let iter = 10 in
let outc = open_out "export.txt" in
for k = 1 to iter do
	plot_coll k outc  @@
	main n m
done;

Printf.fprintf outc "sh = max(-1,len(a1)";

for k = 2 to iter do
	Printf.fprintf outc "%s" (",len(a"^(string_of_int k)^")");
done;

Printf.fprintf outc ")\ndef zerone(x):\n\tif x < 1:\n\t\treturn 1\n\treturn x\na = np.zeros(sh, dtype=np.float64)\n";

for k = 2 to iter do
    Printf.fprintf outc "%s" ("a"^(string_of_int k)^".resize(sh, refcheck=False)\na"^(string_of_int k)^" = np.vectorize(zerone)(a"^(string_of_int k)^")\na += a"^(string_of_int k)^"\n") 
done;

Printf.fprintf outc "%s" ("\na = a / "^(string_of_int iter)^"\n#eV(a[0],sh)\n");
Printf.fprintf outc "\ngraph(a)"



