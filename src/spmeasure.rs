
use std::collections::HashMap;
use std::f64;

use crate::error;
use crate::matrices;

/* Substitution scoring matrix */
static mut M : Vec<f64> = Vec::new();

pub fn sum_of_pairs( site_list : &Vec<String>, /*weight_list : &Vec<f64>,*/ arg_m : &String ) -> Vec<f64>
{
	let num_site : usize = ( *site_list ).len();

	/* 
	 * Amino acid index for picking the substitution scoring matrix elements.
	 * aa_index = Amino acid index
	 * aa_list  = Amino acid order of substitution scoring matrices
	 */
	let mut aa_index : HashMap<char, usize> = HashMap::new();
	let     aa_list  : Vec<char> = "ARNDCQEGHILKMFPSTWYV".chars().collect();
	for i in 0 .. 20 {
		aa_index.insert( aa_list[ i ], i );
	}
	//println!( "aa_index : {:?}", aa_index );

	/* Make a substitution scoring matrix. */
	unsafe { 
		/* Define a  substitution scoring matrix. */
		M = matrices::define_matrix( arg_m );
		M.shrink_to_fit();
		//println!( "{:?}", M );

		/* Check the amino acid index in the substitution scoring matrix. */
		/*
		for i in aa_list.iter() {
			for j in aa_list.iter() {
				println!( "M[ {}, {} ] = {}", *i, *j, M[ ( aa_index )[ i ] * 20 + ( aa_index )[ j ] ] );
			}
		}
		*/

		/* Normalise the scoring matrix based on Karlin-like method */
		M = normalize_matrix( /* &aa_list, &aa_index */ );
		/*
		for a in aa_list.iter() {
			for b in aa_list.iter() {
				println!( "Normalized M[ {}, {} ] = {:.3}", *a, *b, M[ aa_index[ a ] * 20 + aa_index[ b ] ] );
			}
		}
		*/
		println!( "\nMaxima of the matrix : {:.3}", M.iter().fold( 0.0 / 0.0, | m, v | v.max( m ) ) );
		println!( "Minima of the matrix : {:.3}\n", M.iter().fold( 0.0 / 0.0, | m, v | v.min( m ) ) );

		/* Check whether the matrix is diagonal. */
		check_matrix_diag();

	}

	let mut sp_score_list : Vec<f64> = vec![ 0.0; num_site ];
	for i in 0 .. num_site {
		let sp_score : f64 = calc_sp( &( *site_list )[ i ] /*, weight_list*/, &aa_index );
		//println!( "SP score of site {} : {:.3}", i, sp_score );
		sp_score_list[ i ] += sp_score;
	}

	sp_score_list
}

fn normalize_matrix( /* aa_list : &Vec<char>, aa_index : &HashMap<char, usize> */ ) -> Vec<f64>
{
	let mut norm_matrix : Vec<f64> = vec![ 0.0; 400 ];

	/*
	 * Normalise the scoring matrix based on Karlin-like method.
	 * It ensures that M(a,a) is 1.
	 * max = Maxima of the matrix
	 * min = Minima of the matrix
	 */
	unsafe {
		//let mmax : f64 = M.iter().fold( 0.0 / 0.0, | m, v | v.max( m ) );
		//let mmin : f64 = M.iter().fold( 0.0 / 0.0, | m, v | v.min( m ) );
		//println!( "{}, {}", mmax, mmin );
		for a in 0 .. 20 {
			for b in 0 .. 20 {
				let mat_aa : f64 = M[ a * 20 + a ];
				let mat_bb : f64 = M[ b * 20 + b ];
				let mat_ab : f64 = M[ a * 20 + b ];
				norm_matrix[ a * 20 + b ] = mat_ab / ( mat_aa * mat_bb ).sqrt();
			}
		}
	}

	/* Check the size of the matrix. */
	if norm_matrix.len() != 400 {
		error::error_bomb( "mat_not_20*20" );
	}

	norm_matrix
}

fn calc_sp( site : &String /*,weight_list : &Vec<f64>*/, aa_index : &HashMap<char, usize> ) -> f64
{
	let     char_list : Vec<char> = ( *site ).chars().collect();
	let mut sp_score  : f64 = 0.0;
	let     site_len  : usize = char_list.len();

	unsafe {
		/* The minimum element of the normalized matrix is given as gap penalty. */
		let mat_min : f64 = M.iter().fold( 0.0 / 0.0, | m, v | v.min( m ) );

		/*
		 * Calculate residue conservation using Sum-of-pairs measure.
		 * If ( i = gap ) or ( j = gap ), the minimum element of the matrix is given as gap penalty.
		 * i          = Site i
		 * j          = Site j
		 * site_len   = The length of a site
		 * char_list  = Vec<char> of the site (&String)
		 * mat_min    = The minima of the normalized scoring matrix
		 * aa_index[] = The order of amino acids
		 * M[]        = The normalized scoring matrix
		 * sp_score   = Conservation score
		 */
		for i in 0 .. ( site_len - 1 ) {
			for j in ( i + 1 ) .. site_len {
				if char_list[ i ] == '-' {
					//println!( "M[ -, {} ] = {:.3}", char_list[ j ], mat_min );
					sp_score += mat_min;
				} else if char_list[ j ] == '-' {
					//println!( "M[ {}, - ] = {:.3}", char_list[ i ], mat_min );
					sp_score += mat_min;
				} else {
					//println!( "M[ {}, {} ] = {:.3}", char_list[ i ], char_list[ j ], M[ ( *aa_index )[ &char_list[ i ] ] * 20 + ( *aa_index )[ &char_list[ j ] ] ] );
					sp_score += M[ ( *aa_index )[ &char_list[ i ] ] * 20 + ( *aa_index )[ &char_list[ j ] ] ];
				}
			}
		}
	}

	sp_score *= 2.0 / ( site_len * ( site_len - 1 ) ) as f64;

	sp_score
}

fn check_matrix_diag()
{
	unsafe {
		for a in 0 .. 20 {
			for b in 0 .. 20 {
				if M[ a * 20 + b ] != M[ b * 20 + a] {
					error::error_bomb( "mat_not_diag" );
				}
			}
		}

	}

}
