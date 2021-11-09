
pub fn define_matrix( arg_m : &String ) -> Vec<f64>
{
	/* 20 × 20 elements */
	let mut _matrix : Vec<f64> = vec![ 0.0; 400 ];

	let _blosum45 : Vec<f64> = vec![ // BLOSUM45
	 5.0, -2.0, -1.0, -2.0, -1.0, -1.0, -1.0,  0.0, -2.0, -1.0, -1.0, -1.0, -1.0, -2.0, -1.0,  1.0,  0.0, -2.0, -2.0,  0.0, // A
	-2.0,  7.0,  0.0, -1.0, -3.0,  1.0,  0.0, -2.0,  0.0, -3.0, -2.0,  3.0, -1.0, -2.0, -2.0, -1.0, -1.0, -2.0, -1.0, -2.0, // R
	-1.0,  0.0,  6.0,  2.0, -2.0,  0.0,  0.0,  0.0,  1.0, -2.0, -3.0,  0.0, -2.0, -2.0, -2.0,  1.0,  0.0, -4.0, -2.0, -3.0, // N
	-2.0, -1.0,  2.0,  7.0, -3.0,  0.0,  2.0, -1.0,  0.0, -4.0, -3.0,  0.0, -3.0, -4.0, -1.0,  0.0, -1.0, -4.0, -2.0, -3.0, // D
	-1.0, -3.0, -2.0, -3.0, 12.0, -3.0, -3.0, -3.0, -3.0, -3.0, -2.0, -3.0, -2.0, -2.0, -4.0, -1.0, -1.0, -5.0, -3.0, -1.0, // C
	-1.0,  1.0,  0.0,  0.0, -3.0,  6.0,  2.0, -2.0,  1.0, -2.0, -2.0,  1.0,  0.0, -4.0, -1.0,  0.0, -1.0, -2.0, -1.0, -3.0, // Q
	-1.0,  0.0,  0.0,  2.0, -3.0,  2.0,  6.0, -2.0,  0.0, -3.0, -2.0,  1.0, -2.0, -3.0,  0.0,  0.0, -1.0, -3.0, -2.0, -3.0, // E
	 0.0, -2.0,  0.0, -1.0, -3.0, -2.0, -2.0,  7.0, -2.0, -4.0, -3.0, -2.0, -2.0, -3.0, -2.0,  0.0, -2.0, -2.0, -3.0, -3.0, // G
	-2.0,  0.0,  1.0,  0.0, -3.0,  1.0,  0.0, -2.0, 10.0, -3.0, -2.0, -1.0,  0.0, -2.0, -2.0, -1.0, -2.0, -3.0,  2.0, -3.0, // H
	-1.0, -3.0, -2.0, -4.0, -3.0, -2.0, -3.0, -4.0, -3.0,  5.0,  2.0, -3.0,  2.0,  0.0, -2.0, -2.0, -1.0, -2.0,  0.0,  3.0, // I
	-1.0, -2.0, -3.0, -3.0, -2.0, -2.0, -2.0, -3.0, -2.0,  2.0,  5.0, -3.0,  2.0,  1.0, -3.0, -3.0, -1.0, -2.0,  0.0,  1.0, // L
	-1.0,  3.0,  0.0,  0.0, -3.0,  1.0,  1.0, -2.0, -1.0, -3.0, -3.0,  5.0, -1.0, -3.0, -1.0, -1.0, -1.0, -2.0, -1.0, -2.0, // K
	-1.0, -1.0, -2.0, -3.0, -2.0,  0.0, -2.0, -2.0,  0.0,  2.0,  2.0, -1.0,  6.0,  0.0, -2.0, -2.0, -1.0, -2.0,  0.0,  1.0, // M
	-2.0, -2.0, -2.0, -4.0, -2.0, -4.0, -3.0, -3.0, -2.0,  0.0,  1.0, -3.0,  0.0,  8.0, -3.0, -2.0, -1.0,  1.0,  3.0,  0.0, // F
	-1.0, -2.0, -2.0, -1.0, -4.0, -1.0,  0.0, -2.0, -2.0, -2.0, -3.0, -1.0, -2.0, -3.0,  9.0, -1.0, -1.0, -3.0, -3.0, -3.0, // P
	 1.0, -1.0,  1.0,  0.0, -1.0,  0.0,  0.0,  0.0, -1.0, -2.0, -3.0, -1.0, -2.0, -2.0, -1.0,  4.0,  2.0, -4.0, -2.0, -1.0, // S
	 0.0, -1.0,  0.0, -1.0, -1.0, -1.0, -1.0, -2.0, -2.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0,  2.0,  5.0, -3.0, -1.0,  0.0, // T
	-2.0, -2.0, -4.0, -4.0, -5.0, -2.0, -3.0, -2.0, -3.0, -2.0, -2.0, -2.0, -2.0,  1.0, -3.0, -4.0, -3.0, 15.0,  3.0, -3.0, // W
	-2.0, -1.0, -2.0, -2.0, -3.0, -1.0, -2.0, -3.0,  2.0,  0.0,  0.0, -1.0,  0.0,  3.0, -3.0, -2.0, -1.0,  3.0,  8.0, -1.0, // Y
	 0.0, -2.0, -3.0, -3.0, -1.0, -3.0, -3.0, -3.0, -3.0,  3.0,  1.0, -2.0,  1.0,  0.0, -3.0, -1.0,  0.0, -3.0, -1.0,  5.0  // V
	// A,    R,    N,    D,    C,    Q,    E,    G,    H,    I,    L,    K,    M,    F,    P,    S,    T,    W,    Y,    V
	];

	let _blosum50 : Vec<f64> = vec![ // BLOSUM50
	 5.0, -2.0, -1.0, -2.0, -1.0, -1.0, -1.0,  0.0, -2.0, -1.0, -2.0, -1.0, -1.0, -3.0, -1.0,  1.0,  0.0, -3.0, -2.0,  0.0, // A
	-2.0,  7.0, -1.0, -2.0, -4.0,  1.0,  0.0, -3.0,  0.0, -4.0, -3.0,  3.0, -2.0, -3.0, -3.0, -1.0, -1.0, -3.0, -1.0, -3.0, // R
	-1.0, -1.0,  7.0,  2.0, -2.0,  0.0,  0.0,  0.0,  1.0, -3.0, -4.0,  0.0, -2.0, -4.0, -2.0,  1.0,  0.0, -4.0, -2.0, -3.0, // N
	-2.0, -2.0,  2.0,  8.0, -4.0,  0.0,  2.0, -1.0, -1.0, -4.0, -4.0, -1.0, -4.0, -5.0, -1.0,  0.0, -1.0, -5.0, -3.0, -4.0, // D
	-1.0, -4.0, -2.0, -4.0, 13.0, -3.0, -3.0, -3.0, -3.0, -2.0, -2.0, -3.0, -2.0, -2.0, -4.0, -1.0, -1.0, -5.0, -3.0, -1.0, // C
	-1.0,  1.0,  0.0,  0.0, -3.0,  7.0,  2.0, -2.0,  1.0, -3.0, -2.0,  2.0,  0.0, -4.0, -1.0,  0.0, -1.0, -1.0, -1.0, -3.0, // Q
	-1.0,  0.0,  0.0,  2.0, -3.0,  2.0,  6.0, -3.0,  0.0, -4.0, -3.0,  1.0, -2.0, -3.0, -1.0, -1.0, -1.0, -3.0, -2.0, -3.0, // E
	 0.0, -3.0,  0.0, -1.0, -3.0, -2.0, -3.0,  8.0, -2.0, -4.0, -4.0, -2.0, -3.0, -4.0, -2.0,  0.0, -2.0, -3.0, -3.0, -4.0, // G
	-2.0,  0.0,  1.0, -1.0, -3.0,  1.0,  0.0, -2.0, 10.0, -4.0, -3.0,  0.0, -1.0, -1.0, -2.0, -1.0, -2.0, -3.0,  2.0, -4.0, // H
	-1.0, -4.0, -3.0, -4.0, -2.0, -3.0, -4.0, -4.0, -4.0,  5.0,  2.0, -3.0,  2.0,  0.0, -3.0, -3.0, -1.0, -3.0, -1.0,  4.0, // I
	-2.0, -3.0, -4.0, -4.0, -2.0, -2.0, -3.0, -4.0, -3.0,  2.0,  5.0, -3.0,  3.0,  1.0, -4.0, -3.0, -1.0, -2.0, -1.0,  1.0, // L
	-1.0,  3.0,  0.0, -1.0, -3.0,  2.0,  1.0, -2.0,  0.0, -3.0, -3.0,  6.0, -2.0, -4.0, -1.0,  0.0, -1.0, -3.0, -2.0, -3.0, // K
	-1.0, -2.0, -2.0, -4.0, -2.0,  0.0, -2.0, -3.0, -1.0,  2.0,  3.0, -2.0,  7.0,  0.0, -3.0, -2.0, -1.0, -1.0,  0.0,  1.0, // M
	-3.0, -3.0, -4.0, -5.0, -2.0, -4.0, -3.0, -4.0, -1.0,  0.0,  1.0, -4.0,  0.0,  8.0, -4.0, -3.0, -2.0,  1.0,  4.0, -1.0, // F
	-1.0, -3.0, -2.0, -1.0, -4.0, -1.0, -1.0, -2.0, -2.0, -3.0, -4.0, -1.0, -3.0, -4.0, 10.0, -1.0, -1.0, -4.0, -3.0, -3.0, // P
	 1.0, -1.0,  1.0,  0.0, -1.0,  0.0, -1.0,  0.0, -1.0, -3.0, -3.0,  0.0, -2.0, -3.0, -1.0,  5.0,  2.0, -4.0, -2.0, -2.0, // S
	 0.0, -1.0,  0.0, -1.0, -1.0, -1.0, -1.0, -2.0, -2.0, -1.0, -1.0, -1.0, -1.0, -2.0, -1.0,  2.0,  5.0, -3.0, -2.0,  0.0, // T
	-3.0, -3.0, -4.0, -5.0, -5.0, -1.0, -3.0, -3.0, -3.0, -3.0, -2.0, -3.0, -1.0,  1.0, -4.0, -4.0, -3.0, 15.0,  2.0, -3.0, // W
	-2.0, -1.0, -2.0, -3.0, -3.0, -1.0, -2.0, -3.0,  2.0, -1.0, -1.0, -2.0,  0.0,  4.0, -3.0, -2.0, -2.0,  2.0,  8.0, -1.0, // Y
	 0.0, -3.0, -3.0, -4.0, -1.0, -3.0, -3.0, -4.0, -4.0,  4.0,  1.0, -3.0,  1.0, -1.0, -3.0, -2.0,  0.0, -3.0, -1.0,  5.0  // V
	// A,    R,    N,    D,    C,    Q,    E,    G,    H,    I,    L,    K,    M,    F,    P,    S,    T,    W,    Y,    V
	];

	let _blosum62 : Vec<f64> = vec![ // BLOSUM62
	 4.0, -1.0, -2.0, -2.0,  0.0, -1.0, -1.0,  0.0, -2.0, -1.0, -1.0, -1.0, -1.0, -2.0, -1.0,  1.0,  0.0, -3.0, -2.0,  0.0, // A
	-1.0,  5.0,  0.0, -2.0, -3.0,  1.0,  0.0, -2.0,  0.0, -3.0, -2.0,  2.0, -1.0, -3.0, -2.0, -1.0, -1.0, -3.0, -2.0, -3.0, // R
	-2.0,  0.0,  6.0,  1.0, -3.0,  0.0,  0.0,  0.0,  1.0, -3.0, -3.0,  0.0, -2.0, -3.0, -2.0,  1.0,  0.0, -4.0, -2.0, -3.0, // N
	-2.0, -2.0,  1.0,  6.0, -3.0,  0.0,  2.0, -1.0, -1.0, -3.0, -4.0, -1.0, -3.0, -3.0, -1.0,  0.0, -1.0, -4.0, -3.0, -3.0, // D
	 0.0, -3.0, -3.0, -3.0,  9.0, -3.0, -4.0, -3.0, -3.0, -1.0, -1.0, -3.0, -1.0, -2.0, -3.0, -1.0, -1.0, -2.0, -2.0, -1.0, // C
	-1.0,  1.0,  0.0,  0.0, -3.0,  5.0,  2.0, -2.0,  0.0, -3.0, -2.0,  1.0,  0.0, -3.0, -1.0,  0.0, -1.0, -2.0, -1.0, -2.0, // Q
	-1.0,  0.0,  0.0,  2.0, -4.0,  2.0,  5.0, -2.0,  0.0, -3.0, -3.0,  1.0, -2.0, -3.0, -1.0,  0.0, -1.0, -3.0, -2.0, -2.0, // E
	 0.0, -2.0,  0.0, -1.0, -3.0, -2.0, -2.0,  6.0, -2.0, -4.0, -4.0, -2.0, -3.0, -3.0, -2.0,  0.0, -2.0, -2.0, -3.0, -3.0, // G
	-2.0,  0.0,  1.0, -1.0, -3.0,  0.0,  0.0, -2.0,  8.0, -3.0, -3.0, -1.0, -2.0, -1.0, -2.0, -1.0, -2.0, -2.0,  2.0, -3.0, // H
	-1.0, -3.0, -3.0, -3.0, -1.0, -3.0, -3.0, -4.0, -3.0,  4.0,  2.0, -3.0,  1.0,  0.0, -3.0, -2.0, -1.0, -3.0, -1.0,  3.0, // I
	-1.0, -2.0, -3.0, -4.0, -1.0, -2.0, -3.0, -4.0, -3.0,  2.0,  4.0, -2.0,  2.0,  0.0, -3.0, -2.0, -1.0, -2.0, -1.0,  1.0, // L
	-1.0,  2.0,  0.0, -1.0, -3.0,  1.0,  1.0, -2.0, -1.0, -3.0, -2.0,  5.0, -1.0, -3.0, -1.0,  0.0, -1.0, -3.0, -2.0, -2.0, // K
	-1.0, -1.0, -2.0, -3.0, -1.0,  0.0, -2.0, -3.0, -2.0,  1.0,  2.0, -1.0,  5.0,  0.0, -2.0, -1.0, -1.0, -1.0, -1.0,  1.0, // M
	-2.0, -3.0, -3.0, -3.0, -2.0, -3.0, -3.0, -3.0, -1.0,  0.0,  0.0, -3.0,  0.0,  6.0, -4.0, -2.0, -2.0,  1.0,  3.0, -1.0, // F
	-1.0, -2.0, -2.0, -1.0, -3.0, -1.0, -1.0, -2.0, -2.0, -3.0, -3.0, -1.0, -2.0, -4.0,  7.0, -1.0, -1.0, -4.0, -3.0, -2.0, // P
	 1.0, -1.0,  1.0,  0.0, -1.0,  0.0,  0.0,  0.0, -1.0, -2.0, -2.0,  0.0, -1.0, -2.0, -1.0,  4.0,  1.0, -3.0, -2.0, -2.0, // S
	 0.0, -1.0,  0.0, -1.0, -1.0, -1.0, -1.0, -2.0, -2.0, -1.0, -1.0, -1.0, -1.0, -2.0, -1.0,  1.0,  5.0, -2.0, -2.0,  0.0, // T
	-3.0, -3.0, -4.0, -4.0, -2.0, -2.0, -3.0, -2.0, -2.0, -3.0, -2.0, -3.0, -1.0,  1.0, -4.0, -3.0, -2.0, 11.0,  2.0, -3.0, // W
	-2.0, -2.0, -2.0, -3.0, -2.0, -1.0, -2.0, -3.0,  2.0, -1.0, -1.0, -2.0, -1.0,  3.0, -3.0, -2.0, -2.0,  2.0,  7.0, -1.0, // Y
	 0.0, -3.0, -3.0, -3.0, -1.0, -2.0, -2.0, -3.0, -3.0,  3.0,  1.0, -2.0,  1.0, -1.0, -2.0, -2.0,  0.0, -3.0, -1.0,  4.0  // V
	// A,    R,    N,    D,    C,    Q,    E,    G,    H,    I,    L,    K,    M,    F,    P,    S,    T,    W,    Y,    V
	];

	let _blosum80 : Vec<f64> = vec![ // BLOSUM80
	 5.0, -2.0, -2.0, -2.0, -1.0, -1.0, -1.0,  0.0, -2.0, -2.0, -2.0, -1.0, -1.0, -3.0, -1.0,  1.0,  0.0, -3.0, -2.0,  0.0, // A
	-2.0,  6.0, -1.0, -2.0, -4.0,  1.0, -1.0, -3.0,  0.0, -3.0, -3.0,  2.0, -2.0, -4.0, -2.0, -1.0, -1.0, -4.0, -3.0, -3.0, // R
	-2.0, -1.0,  6.0,  1.0, -3.0,  0.0, -1.0, -1.0,  0.0, -4.0, -4.0,  0.0, -3.0, -4.0, -3.0,  0.0,  0.0, -4.0, -3.0, -4.0, // N
	-2.0, -2.0,  1.0,  6.0, -4.0, -1.0,  1.0, -2.0, -2.0, -4.0, -5.0, -1.0, -4.0, -4.0, -2.0, -1.0, -1.0, -6.0, -4.0, -4.0, // D
	-1.0, -4.0, -3.0, -4.0,  9.0, -4.0, -5.0, -4.0, -4.0, -2.0, -2.0, -4.0, -2.0, -3.0, -4.0, -2.0, -1.0, -3.0, -3.0, -1.0, // C
	-1.0,  1.0,  0.0, -1.0, -4.0,  6.0,  2.0, -2.0,  1.0, -3.0, -3.0,  1.0,  0.0, -4.0, -2.0,  0.0, -1.0, -3.0, -2.0, -3.0, // Q
	-1.0, -1.0, -1.0,  1.0, -5.0,  2.0,  6.0, -3.0,  0.0, -4.0, -4.0,  1.0, -2.0, -4.0, -2.0,  0.0, -1.0, -4.0, -3.0, -3.0, // E
	 0.0, -3.0, -1.0, -2.0, -4.0, -2.0, -3.0,  6.0, -3.0, -5.0, -4.0, -2.0, -4.0, -4.0, -3.0, -1.0, -2.0, -4.0, -4.0, -4.0, // G
	-2.0,  0.0,  0.0, -2.0, -4.0,  1.0,  0.0, -3.0,  8.0, -4.0, -3.0, -1.0, -2.0, -2.0, -3.0, -1.0, -2.0, -3.0,  2.0, -4.0, // H
	-2.0, -3.0, -4.0, -4.0, -2.0, -3.0, -4.0, -5.0, -4.0,  5.0,  1.0, -3.0,  1.0, -1.0, -4.0, -3.0, -1.0, -3.0, -2.0,  3.0, // I
	-2.0, -3.0, -4.0, -5.0, -2.0, -3.0, -4.0, -4.0, -3.0,  1.0,  4.0, -3.0,  2.0,  0.0, -3.0, -3.0, -2.0, -2.0, -2.0,  1.0, // L
	-1.0,  2.0,  0.0, -1.0, -4.0,  1.0,  1.0, -2.0, -1.0, -3.0, -3.0,  5.0, -2.0, -4.0, -1.0, -1.0, -1.0, -4.0, -3.0, -3.0, // K
	-1.0, -2.0, -3.0, -4.0, -2.0,  0.0, -2.0, -4.0, -2.0,  1.0,  2.0, -2.0,  6.0,  0.0, -3.0, -2.0, -1.0, -2.0, -2.0,  1.0, // M
	-3.0, -4.0, -4.0, -4.0, -3.0, -4.0, -4.0, -4.0, -2.0, -1.0,  0.0, -4.0,  0.0,  6.0, -4.0, -3.0, -2.0,  0.0,  3.0, -1.0, // F
	-1.0, -2.0, -3.0, -2.0, -4.0, -2.0, -2.0, -3.0, -3.0, -4.0, -3.0, -1.0, -3.0, -4.0,  8.0, -1.0, -2.0, -5.0, -4.0, -3.0, // P
	 1.0, -1.0,  0.0, -1.0, -2.0,  0.0,  0.0, -1.0, -1.0, -3.0, -3.0, -1.0, -2.0, -3.0, -1.0,  5.0,  1.0, -4.0, -2.0, -2.0, // S
	 0.0, -1.0,  0.0, -1.0, -1.0, -1.0, -1.0, -2.0, -2.0, -1.0, -2.0, -1.0, -1.0, -2.0, -2.0,  1.0,  5.0, -4.0, -2.0,  0.0, // T
	-3.0, -4.0, -4.0, -6.0, -3.0, -3.0, -4.0, -4.0, -3.0, -3.0, -2.0, -4.0, -2.0,  0.0, -5.0, -4.0, -4.0, 11.0,  2.0, -3.0, // W
	-2.0, -3.0, -3.0, -4.0, -3.0, -2.0, -3.0, -4.0,  2.0, -2.0, -2.0, -3.0, -2.0,  3.0, -4.0, -2.0, -2.0,  2.0,  7.0, -2.0, // Y
	 0.0, -3.0, -4.0, -4.0, -1.0, -3.0, -3.0, -4.0, -4.0,  3.0,  1.0, -3.0,  1.0, -1.0, -3.0, -2.0,  0.0, -3.0, -2.0,  4.0  // V
	// A,    R,    N,    D,    C,    Q,    E,    G,    H,    I,    L,    K,    M,    F,    P,    S,    T,    W,    Y,    V
	];

	let _blosum90 : Vec<f64> = vec![ // BLOSUM90
	 5.0, -2.0, -2.0, -3.0, -1.0, -1.0, -1.0,  0.0, -2.0, -2.0, -2.0, -1.0, -2.0, -3.0, -1.0,  1.0,  0.0, -4.0, -3.0, -1.0, // A
	-2.0,  6.0, -1.0, -3.0, -5.0,  1.0, -1.0, -3.0,  0.0, -4.0, -3.0,  2.0, -2.0, -4.0, -3.0, -1.0, -2.0, -4.0, -3.0, -3.0, // R
	-2.0, -1.0,  7.0,  1.0, -4.0,  0.0, -1.0, -1.0,  0.0, -4.0, -4.0,  0.0, -3.0, -4.0, -3.0,  0.0,  0.0, -5.0, -3.0, -4.0, // N
	-3.0, -3.0,  1.0,  7.0, -5.0, -1.0,  1.0, -2.0, -2.0, -5.0, -5.0, -1.0, -4.0, -5.0, -3.0, -1.0, -2.0, -6.0, -4.0, -5.0, // D
	-1.0, -5.0, -4.0, -5.0,  9.0, -4.0, -6.0, -4.0, -5.0, -2.0, -2.0, -4.0, -2.0, -3.0, -4.0, -2.0, -2.0, -4.0, -4.0, -2.0, // C
	-1.0,  1.0,  0.0, -1.0, -4.0,  7.0,  2.0, -3.0,  1.0, -4.0, -3.0,  1.0,  0.0, -4.0, -2.0, -1.0, -1.0, -3.0, -3.0, -3.0, // Q
	-1.0, -1.0, -1.0,  1.0, -6.0,  2.0,  6.0, -3.0, -1.0, -4.0, -4.0,  0.0, -3.0, -5.0, -2.0, -1.0, -1.0, -5.0, -4.0, -3.0, // E
	 0.0, -3.0, -1.0, -2.0, -4.0, -3.0, -3.0,  6.0, -3.0, -5.0, -5.0, -2.0, -4.0, -5.0, -3.0, -1.0, -3.0, -4.0, -5.0, -5.0, // G
	-2.0,  0.0,  0.0, -2.0, -5.0,  1.0, -1.0, -3.0,  8.0, -4.0, -4.0, -1.0, -3.0, -2.0, -3.0, -2.0, -2.0, -3.0,  1.0, -4.0, // H
	-2.0, -4.0, -4.0, -5.0, -2.0, -4.0, -4.0, -5.0, -4.0,  5.0,  1.0, -4.0,  1.0, -1.0, -4.0, -3.0, -1.0, -4.0, -2.0,  3.0, // I
	-2.0, -3.0, -4.0, -5.0, -2.0, -3.0, -4.0, -5.0, -4.0,  1.0,  5.0, -3.0,  2.0,  0.0, -4.0, -3.0, -2.0, -3.0, -2.0,  0.0, // L
	-1.0,  2.0,  0.0, -1.0, -4.0,  1.0,  0.0, -2.0, -1.0, -4.0, -3.0,  6.0, -2.0, -4.0, -2.0, -1.0, -1.0, -5.0, -3.0, -3.0, // K
	-2.0, -2.0, -3.0, -4.0, -2.0,  0.0, -3.0, -4.0, -3.0,  1.0,  2.0, -2.0,  7.0, -1.0, -3.0, -2.0, -1.0, -2.0, -2.0,  0.0, // M
	-3.0, -4.0, -4.0, -5.0, -3.0, -4.0, -5.0, -5.0, -2.0, -1.0,  0.0, -4.0, -1.0,  7.0, -4.0, -3.0, -3.0,  0.0,  3.0, -2.0, // F
	-1.0, -3.0, -3.0, -3.0, -4.0, -2.0, -2.0, -3.0, -3.0, -4.0, -4.0, -2.0, -3.0, -4.0,  8.0, -2.0, -2.0, -5.0, -4.0, -3.0, // P
	 1.0, -1.0,  0.0, -1.0, -2.0, -1.0, -1.0, -1.0, -2.0, -3.0, -3.0, -1.0, -2.0, -3.0, -2.0,  5.0,  1.0, -4.0, -3.0, -2.0, // S
	 0.0, -2.0,  0.0, -2.0, -2.0, -1.0, -1.0, -3.0, -2.0, -1.0, -2.0, -1.0, -1.0, -3.0, -2.0,  1.0,  6.0, -4.0, -2.0, -1.0, // T
	-4.0, -4.0, -5.0, -6.0, -4.0, -3.0, -5.0, -4.0, -3.0, -4.0, -3.0, -5.0, -2.0,  0.0, -5.0, -4.0, -4.0, 11.0,  2.0, -3.0, // W
	-3.0, -3.0, -3.0, -4.0, -4.0, -3.0, -4.0, -5.0,  1.0, -2.0, -2.0, -3.0, -2.0,  3.0, -4.0, -3.0, -2.0,  2.0,  8.0, -3.0, // Y
	-1.0, -3.0, -4.0, -5.0, -2.0, -3.0, -3.0, -5.0, -4.0,  3.0,  0.0, -3.0,  0.0, -2.0, -3.0, -2.0, -1.0, -3.0, -3.0,  5.0  // V
	// A,    R,    N,    D,    C,    Q,    E,    G,    H,    I,    L,    K,    M,    F,    P,    S,    T,    W,    Y,    V
	];

	let _pam30 : Vec<f64> = vec![ // PAM30
	  6.0,  -7.0,  -4.0,  -3.0,  -6.0,  -4.0,  -2.0,  -2.0,  -7.0,  -5.0,  -6.0,  -7.0,  -5.0,  -8.0,  -2.0,   0.0,  -1.0, -13.0,  -8.0,  -2.0, // A
	 -7.0,   8.0,  -6.0, -10.0,  -8.0,  -2.0,  -9.0,  -9.0,  -2.0,  -5.0,  -8.0,   0.0,  -4.0,  -9.0,  -4.0,  -3.0,  -6.0,  -2.0, -10.0,  -8.0, // R
	 -4.0,  -6.0,   8.0,   2.0, -11.0,  -3.0,  -2.0,  -3.0,   0.0,  -5.0,  -7.0,  -1.0,  -9.0,  -9.0,  -6.0,   0.0,  -2.0,  -8.0,  -4.0,  -8.0, // N
	 -3.0, -10.0,   2.0,   8.0, -14.0,  -2.0,   2.0,  -3.0,  -4.0,  -7.0, -12.0,  -4.0, -11.0, -15.0,  -8.0,  -4.0,  -5.0, -15.0, -11.0,  -8.0, // D
	 -6.0,  -8.0, -11.0, -14.0,  10.0, -14.0, -14.0,  -9.0,  -7.0,  -6.0, -15.0, -14.0, -13.0, -13.0,  -8.0,  -3.0,  -8.0, -15.0,  -4.0,  -6.0, // C
	 -4.0,  -2.0,  -3.0,  -2.0, -14.0,   8.0,   1.0,  -7.0,   1.0,  -8.0,  -5.0,  -3.0,  -4.0, -13.0,  -3.0,  -5.0,  -5.0, -13.0, -12.0,  -7.0, // Q
	 -2.0,  -9.0,  -2.0,   2.0, -14.0,   1.0,   8.0,  -4.0,  -5.0,  -5.0,  -9.0,  -4.0,  -7.0, -14.0,  -5.0,  -4.0,  -6.0, -17.0,  -8.0,  -6.0, // E
	 -2.0,  -9.0,  -3.0,  -3.0,  -9.0,  -7.0,  -4.0,   6.0,  -9.0, -11.0, -10.0,  -7.0,  -8.0,  -9.0,  -6.0,  -2.0,  -6.0, -15.0, -14.0,  -5.0, // G
	 -7.0,  -2.0,   0.0,  -4.0,  -7.0,   1.0,  -5.0,  -9.0,   9.0,  -9.0,  -6.0,  -6.0, -10.0,  -6.0,  -4.0,  -6.0,  -7.0,  -7.0,  -3.0,  -6.0, // H
	 -5.0,  -5.0,  -5.0,  -7.0,  -6.0,  -8.0,  -5.0, -11.0,  -9.0,   8.0,  -1.0,  -6.0,  -1.0,  -2.0,  -8.0,  -7.0,  -2.0, -14.0,  -6.0,   2.0, // I
	 -6.0,  -8.0,  -7.0, -12.0, -15.0,  -5.0,  -9.0, -10.0,  -6.0,  -1.0,   7.0,  -8.0,   1.0,  -3.0,  -7.0,  -8.0,  -7.0,  -6.0,  -7.0,  -2.0, // L
	 -7.0,   0.0,  -1.0,  -4.0, -14.0,  -3.0,  -4.0,  -7.0,  -6.0,  -6.0,  -8.0,   7.0,  -2.0, -14.0,  -6.0,  -4.0,  -3.0, -12.0,  -9.0,  -9.0, // K
	 -5.0,  -4.0,  -9.0, -11.0, -13.0,  -4.0,  -7.0,  -8.0, -10.0,  -1.0,   1.0,  -2.0,  11.0,  -4.0,  -8.0,  -5.0,  -4.0, -13.0, -11.0,  -1.0, // M
	 -8.0,  -9.0,  -9.0, -15.0, -13.0, -13.0, -14.0,  -9.0,  -6.0,  -2.0,  -3.0, -14.0,  -4.0,   9.0, -10.0,  -6.0,  -9.0,  -4.0,   2.0,  -8.0, // F
	 -2.0,  -4.0,  -6.0,  -8.0,  -8.0,  -3.0,  -5.0,  -6.0,  -4.0,  -8.0,  -7.0,  -6.0,  -8.0, -10.0,   8.0,  -2.0,  -4.0, -14.0, -13.0,  -6.0, // P
	  0.0,  -3.0,   0.0,  -4.0,  -3.0,  -5.0,  -4.0,  -2.0,  -6.0,  -7.0,  -8.0,  -4.0,  -5.0,  -6.0,  -2.0,   6.0,   0.0,  -5.0,  -7.0,  -6.0, // S
	 -1.0,  -6.0,  -2.0,  -5.0,  -8.0,  -5.0,  -6.0,  -6.0,  -7.0,  -2.0,  -7.0,  -3.0,  -4.0,  -9.0,  -4.0,   0.0,   7.0, -13.0,  -6.0,  -3.0, // T
	-13.0,  -2.0,  -8.0, -15.0, -15.0, -13.0, -17.0, -15.0,  -7.0, -14.0,  -6.0, -12.0, -13.0,  -4.0, -14.0,  -5.0, -13.0,  13.0,  -5.0, -15.0, // W
	 -8.0, -10.0,  -4.0, -11.0,  -4.0, -12.0,  -8.0, -14.0,  -3.0,  -6.0,  -7.0,  -9.0, -11.0,   2.0, -13.0,  -7.0,  -6.0,  -5.0,  10.0,  -7.0, // Y
	 -2.0,  -8.0,  -8.0,  -8.0,  -6.0,  -7.0,  -6.0,  -5.0,  -6.0,   2.0,  -2.0,  -9.0,  -1.0,  -8.0,  -6.0,  -6.0,  -3.0, -15.0,  -7.0,   7.0  // V
	//  A,     R,     N,     D,     C,     Q,     E,     G,     H,     I,     L,     K,     M,     F,     P,     S,     T,     W,     Y,     V
	];

	let _pam70 : Vec<f64> = vec![ // PAM70
	 5.0,  -4.0,  -2.0,  -1.0,  -4.0,  -2.0,  -1.0,   0.0,  -4.0,  -2.0,  -4.0,  -4.0,  -3.0,  -6.0,   0.0,   1.0,   1.0,  -9.0,  -5.0,  -1.0, // A
	-4.0,   8.0,  -3.0,  -6.0,  -5.0,   0.0,  -5.0,  -6.0,   0.0,  -3.0,  -6.0,   2.0,  -2.0,  -7.0,  -2.0,  -1.0,  -4.0,   0.0,  -7.0,  -5.0, // R
	-2.0,  -3.0,   6.0,   3.0,  -7.0,  -1.0,   0.0,  -1.0,   1.0,  -3.0,  -5.0,   0.0,  -5.0,  -6.0,  -3.0,   1.0,   0.0,  -6.0,  -3.0,  -5.0, // N
	-1.0,  -6.0,   3.0,   6.0,  -9.0,   0.0,   3.0,  -1.0,  -1.0,  -5.0,  -8.0,  -2.0,  -7.0, -10.0,  -4.0,  -1.0,  -2.0, -10.0,  -7.0,  -5.0, // D
	-4.0,  -5.0,  -7.0,  -9.0,   9.0,  -9.0,  -9.0,  -6.0,  -5.0,  -4.0, -10.0,  -9.0,  -9.0,  -8.0,  -5.0,  -1.0,  -5.0, -11.0,  -2.0,  -4.0, // C
	-2.0,   0.0,  -1.0,   0.0,  -9.0,   7.0,   2.0,  -4.0,   2.0,  -5.0,  -3.0,  -1.0,  -2.0,  -9.0,  -1.0,  -3.0,  -3.0,  -8.0,  -8.0,  -4.0, // Q
	-1.0,  -5.0,   0.0,   3.0,  -9.0,   2.0,   6.0,  -2.0,  -2.0,  -4.0,  -6.0,  -2.0,  -4.0,  -9.0,  -3.0,  -2.0,  -3.0, -11.0,  -6.0,  -4.0, // E
	 0.0,  -6.0,  -1.0,  -1.0,  -6.0,  -4.0,  -2.0,   6.0,  -6.0,  -6.0,  -7.0,  -5.0,  -6.0,  -7.0,  -3.0,   0.0,  -3.0, -10.0,  -9.0,  -3.0, // G
	-4.0,   0.0,   1.0,  -1.0,  -5.0,   2.0,  -2.0,  -6.0,   8.0,  -6.0,  -4.0,  -3.0,  -6.0,  -4.0,  -2.0,  -3.0,  -4.0,  -5.0,  -1.0,  -4.0, // H
	-2.0,  -3.0,  -3.0,  -5.0,  -4.0,  -5.0,  -4.0,  -6.0,  -6.0,   7.0,   1.0,  -4.0,   1.0,   0.0,  -5.0,  -4.0,  -1.0,  -9.0,  -4.0,   3.0, // I
	-4.0,  -6.0,  -5.0,  -8.0, -10.0,  -3.0,  -6.0,  -7.0,  -4.0,   1.0,   6.0,  -5.0,   2.0,  -1.0,  -5.0,  -6.0,  -4.0,  -4.0,  -4.0,   0.0, // L
	-4.0,   2.0,   0.0,  -2.0,  -9.0,  -1.0,  -2.0,  -5.0,  -3.0,  -4.0,  -5.0,   6.0,   0.0,  -9.0,  -4.0,  -2.0,  -1.0,  -7.0,  -7.0,  -6.0, // K
	-3.0,  -2.0,  -5.0,  -7.0,  -9.0,  -2.0,  -4.0,  -6.0,  -6.0,   1.0,   2.0,   0.0,  10.0,  -2.0,  -5.0,  -3.0,  -2.0,  -8.0,  -7.0,   0.0, // M
	-6.0,  -7.0,  -6.0, -10.0,  -8.0,  -9.0,  -9.0,  -7.0,  -4.0,   0.0,  -1.0,  -9.0,  -2.0,   8.0,  -7.0,  -4.0,  -6.0,  -2.0,   4.0,  -5.0, // F
	 0.0,  -2.0,  -3.0,  -4.0,  -5.0,  -1.0,  -3.0,  -3.0,  -2.0,  -5.0,  -5.0,  -4.0,  -5.0,  -7.0,   7.0,   0.0,  -2.0,  -9.0,  -9.0,  -3.0, // P
	 1.0,  -1.0,   1.0,  -1.0,  -1.0,  -3.0,  -2.0,   0.0,  -3.0,  -4.0,  -6.0,  -2.0,  -3.0,  -4.0,   0.0,   5.0,   2.0,  -3.0,  -5.0,  -3.0, // S
	 1.0,  -4.0,   0.0,  -2.0,  -5.0,  -3.0,  -3.0,  -3.0,  -4.0,  -1.0,  -4.0,  -1.0,  -2.0,  -6.0,  -2.0,   2.0,   6.0,  -8.0,  -4.0,  -1.0, // T
	-9.0,   0.0,  -6.0, -10.0, -11.0,  -8.0, -11.0, -10.0,  -5.0,  -9.0,  -4.0,  -7.0,  -8.0,  -2.0,  -9.0,  -3.0,  -8.0,  13.0,  -3.0, -10.0, // W
	-5.0,  -7.0,  -3.0,  -7.0,  -2.0,  -8.0,  -6.0,  -9.0,  -1.0,  -4.0,  -4.0,  -7.0,  -7.0,   4.0,  -9.0,  -5.0,  -4.0,  -3.0,   9.0,  -5.0, // Y
	-1.0,  -5.0,  -5.0,  -5.0,  -4.0,  -4.0,  -4.0,  -3.0,  -4.0,   3.0,   0.0,  -6.0,   0.0,  -5.0,  -3.0,  -3.0,  -1.0, -10.0,  -5.0,   6.0  // V
	// A,     R,     N,     D,     C,     Q,     E,     G,     H,     I,     L,     K,     M,     F,     P,     S,     T,     W,     Y,     V
	];

	let _pam250 : Vec<f64> = vec![ // PAM250
	 2.0, -2.0,  0.0,  0.0, -2.0,  0.0,  0.0,  1.0, -1.0, -1.0, -2.0, -1.0, -1.0, -3.0,  1.0,  1.0,  1.0, -6.0, -3.0,  0.0, // A
	-2.0,  6.0,  0.0, -1.0, -4.0,  1.0, -1.0, -3.0,  2.0, -2.0, -3.0,  3.0,  0.0, -4.0,  0.0,  0.0, -1.0,  2.0, -4.0, -2.0, // R
	 0.0,  0.0,  2.0,  2.0, -4.0,  1.0,  1.0,  0.0,  2.0, -2.0, -3.0,  1.0, -2.0, -3.0,  0.0,  1.0,  0.0, -4.0, -2.0, -2.0, // N
	 0.0, -1.0,  2.0,  4.0, -5.0,  2.0,  3.0,  1.0,  1.0, -2.0, -4.0,  0.0, -3.0, -6.0, -1.0,  0.0,  0.0, -7.0, -4.0, -2.0, // D
	-2.0, -4.0, -4.0, -5.0, 12.0, -5.0, -5.0, -3.0, -3.0, -2.0, -6.0, -5.0, -5.0, -4.0, -3.0,  0.0, -2.0, -8.0,  0.0, -2.0, // C
	 0.0,  1.0,  1.0,  2.0, -5.0,  4.0,  2.0, -1.0,  3.0, -2.0, -2.0,  1.0, -1.0, -5.0,  0.0, -1.0, -1.0, -5.0, -4.0, -2.0, // Q
	 0.0, -1.0,  1.0,  3.0, -5.0,  2.0,  4.0,  0.0,  1.0, -2.0, -3.0,  0.0, -2.0, -5.0, -1.0,  0.0,  0.0, -7.0, -4.0, -2.0, // E
	 1.0, -3.0,  0.0,  1.0, -3.0, -1.0,  0.0,  5.0, -2.0, -3.0, -4.0, -2.0, -3.0, -5.0,  0.0,  1.0,  0.0, -7.0, -5.0, -1.0, // G
	-1.0,  2.0,  2.0,  1.0, -3.0,  3.0,  1.0, -2.0,  6.0, -2.0, -2.0,  0.0, -2.0, -2.0,  0.0, -1.0, -1.0, -3.0,  0.0, -2.0, // H
	-1.0, -2.0, -2.0, -2.0, -2.0, -2.0, -2.0, -3.0, -2.0,  5.0,  2.0, -2.0,  2.0,  1.0, -2.0, -1.0,  0.0, -5.0, -1.0,  4.0, // I
	-2.0, -3.0, -3.0, -4.0, -6.0, -2.0, -3.0, -4.0, -2.0,  2.0,  6.0, -3.0,  4.0,  2.0, -3.0, -3.0, -2.0, -2.0, -1.0,  2.0, // L
	-1.0,  3.0,  1.0,  0.0, -5.0,  1.0,  0.0, -2.0,  0.0, -2.0, -3.0,  5.0,  0.0, -5.0, -1.0,  0.0,  0.0, -3.0, -4.0, -2.0, // K
	-1.0,  0.0, -2.0, -3.0, -5.0, -1.0, -2.0, -3.0, -2.0,  2.0,  4.0,  0.0,  6.0,  0.0, -2.0, -2.0, -1.0, -4.0, -2.0,  2.0, // M
	-3.0, -4.0, -3.0, -6.0, -4.0, -5.0, -5.0, -5.0, -2.0,  1.0,  2.0, -5.0,  0.0,  9.0, -5.0, -3.0, -3.0,  0.0,  7.0, -1.0, // F
	 1.0,  0.0,  0.0, -1.0, -3.0,  0.0, -1.0,  0.0,  0.0, -2.0, -3.0, -1.0, -2.0, -5.0,  6.0,  1.0,  0.0, -6.0, -5.0, -1.0, // P
	 1.0,  0.0,  1.0,  0.0,  0.0, -1.0,  0.0,  1.0, -1.0, -1.0, -3.0,  0.0, -2.0, -3.0,  1.0,  2.0,  1.0, -2.0, -3.0, -1.0, // S
	 1.0, -1.0,  0.0,  0.0, -2.0, -1.0,  0.0,  0.0, -1.0,  0.0, -2.0,  0.0, -1.0, -3.0,  0.0,  1.0,  3.0, -5.0, -3.0,  0.0, // T
	-6.0,  2.0, -4.0, -7.0, -8.0, -5.0, -7.0, -7.0, -3.0, -5.0, -2.0, -3.0, -4.0,  0.0, -6.0, -2.0, -5.0, 17.0,  0.0, -6.0, // W
	-3.0, -4.0, -2.0, -4.0,  0.0, -4.0, -4.0, -5.0,  0.0, -1.0, -1.0, -4.0, -2.0,  7.0, -5.0, -3.0, -3.0,  0.0, 10.0, -2.0, // Y
	 0.0, -2.0, -2.0, -2.0, -2.0, -2.0, -2.0, -1.0, -2.0,  4.0,  2.0, -2.0,  2.0, -1.0, -1.0, -1.0,  0.0, -6.0, -2.0,  4.0  // V
	// A,    R,    N,    D,    C,    Q,    E,    G,    H,    I,    L,    K,    M,    F,    P,    S,    T,    W,    Y,    V
	];

	let _pet91mod : Vec<f64> = vec![ // Modified version of PET91
	15.0, -1.0,  0.0, -1.0, -1.0, -1.0, -1.0,  1.0, -2.0,  0.0, -1.0, -1.0, -1.0, -3.0,  1.0,  1.0,  2.0, -4.0, -3.0,  1.0, // A
	-1.0, 15.0,  0.0, -1.0, -1.0,  2.0,  0.0,  0.0,  2.0, -3.0, -3.0,  4.0, -2.0, -4.0, -1.0, -1.0, -1.0,  0.0, -2.0, -3.0, // R
	 0.0,  0.0, 15.0,  2.0, -1.0,  0.0,  1.0,  0.0,  1.0, -2.0, -3.0,  1.0, -2.0, -3.0, -1.0,  1.0,  1.0, -4.0, -1.0, -2.0, // N
	-1.0, -1.0,  2.0, 15.0, -3.0,  0.0,  4.0,  1.0,  0.0, -3.0, -4.0,  0.0, -3.0, -5.0, -2.0,  0.0, -1.0, -5.0, -2.0, -3.0, // D
	-1.0, -1.0, -1.0, -3.0, 15.0, -3.0, -4.0, -1.0,  0.0, -2.0, -3.0, -3.0, -2.0,  0.0, -2.0,  1.0, -1.0,  1.0,  2.0, -2.0, // C
	-1.0,  2.0,  0.0,  0.0, -3.0, 15.0,  2.0, -1.0,  3.0, -3.0, -2.0,  2.0, -2.0, -4.0,  0.0, -1.0, -1.0, -3.0, -1.0, -3.0, // Q
	-1.0,  0.0,  1.0,  4.0, -4.0,  2.0, 15.0,  1.0,  0.0, -3.0, -4.0,  1.0, -3.0, -5.0, -2.0, -1.0, -1.0, -5.0, -4.0, -2.0, // E
	 1.0,  0.0,  0.0,  1.0, -1.0, -1.0,  1.0, 15.0, -2.0, -3.0, -4.0, -1.0, -3.0, -5.0, -1.0,  1.0,  0.0, -2.0, -4.0, -2.0, // G
	-2.0,  2.0,  1.0,  0.0,  0.0,  3.0,  0.0, -2.0, 15.0, -3.0, -2.0,  1.0, -2.0,  0.0,  0.0, -1.0, -1.0, -3.0,  4.0, -3.0, // H
	 0.0, -3.0, -2.0, -3.0, -2.0, -3.0, -3.0, -3.0, -3.0, 15.0,  2.0, -3.0,  3.0,  0.0, -2.0, -1.0,  1.0, -4.0, -2.0,  4.0, // I
	-1.0, -3.0, -3.0, -4.0, -3.0, -2.0, -4.0, -4.0, -2.0,  2.0, 15.0, -3.0,  3.0,  2.0,  0.0, -2.0, -1.0, -2.0, -1.0,  2.0, // L
	-1.0,  4.0,  1.0,  0.0, -3.0,  2.0,  1.0, -1.0,  1.0, -3.0, -3.0, 15.0, -2.0, -5.0, -2.0, -1.0, -1.0, -3.0, -3.0, -3.0, // K
	-1.0, -2.0, -2.0, -3.0, -2.0, -2.0, -3.0, -3.0, -2.0,  3.0,  3.0, -2.0, 15.0,  0.0, -2.0, -1.0,  0.0, -3.0, -3.0,  2.0, // M
	-3.0, -4.0, -3.0, -5.0,  0.0, -4.0, -5.0, -5.0,  0.0,  0.0,  2.0, -5.0,  0.0, 15.0, -2.0, -2.0, -2.0, -1.0,  5.0,  0.0, // F
	 1.0, -1.0, -1.0, -2.0, -2.0,  0.0, -2.0, -1.0,  0.0, -2.0,  0.0, -2.0, -2.0, -2.0, 15.0,  1.0,  1.0, -5.0, -3.0, -1.0, // P
	 1.0, -1.0,  1.0,  0.0,  1.0, -1.0, -1.0,  1.0, -1.0, -1.0, -2.0, -1.0, -1.0, -2.0,  1.0, 15.0,  1.0, -3.0, -1.0, -1.0, // S
	 2.0, -1.0,  1.0, -1.0, -1.0, -1.0, -1.0,  0.0, -1.0,  1.0, -1.0, -1.0,  0.0, -2.0,  1.0,  1.0, 15.0, -4.0, -3.0,  0.0, // T
	-4.0,  0.0, -4.0, -5.0,  1.0, -3.0, -5.0, -2.0, -3.0, -4.0, -2.0, -3.0, -3.0, -1.0, -5.0, -3.0, -4.0, 15.0,  0.0, -4.0, // W
	-3.0, -2.0, -1.0, -2.0,  2.0, -1.0, -4.0, -4.0,  4.0, -2.0, -1.0, -3.0, -3.0,  5.0, -3.0, -1.0, -3.0,  0.0, 15.0, -3.0, // Y
	 1.0, -3.0, -2.0, -3.0, -2.0, -3.0, -2.0, -2.0, -3.0,  4.0,  2.0, -3.0,  2.0,  0.0, -1.0, -1.0,  0.0, -4.0, -3.0, 15.0  // V
	// A,    R,    N,    D,    C,    Q,    E,    G,    H,    I,    L,    K,    M,    F,    P,    S,    T,    W,    Y,    V
	];

	let _blosum62mod : Vec<f64> = vec![ // Modified version of BLOSUM62
	11.0, -1.0, -2.0, -2.0,  0.0, -1.0, -1.0,  0.0, -2.0, -1.0, -1.0, -1.0, -1.0, -2.0, -1.0,  1.0,  0.0, -3.0, -2.0,  0.0, // A
	-1.0, 11.0,  0.0, -2.0, -3.0,  1.0,  0.0, -2.0,  0.0, -3.0, -2.0,  2.0, -1.0, -3.0, -2.0, -1.0, -1.0, -3.0, -2.0, -3.0, // R
	-2.0,  0.0, 11.0,  1.0, -3.0,  0.0,  0.0,  0.0,  1.0, -3.0, -3.0,  0.0, -2.0, -3.0, -2.0,  1.0,  0.0, -4.0, -2.0, -3.0, // N
	-2.0, -2.0,  1.0, 11.0, -3.0,  0.0,  2.0, -1.0, -1.0, -3.0, -4.0, -1.0, -3.0, -3.0, -1.0,  0.0, -1.0, -4.0, -3.0, -3.0, // D
	 0.0, -3.0, -3.0, -3.0, 11.0, -3.0, -4.0, -3.0, -3.0, -1.0, -1.0, -3.0, -1.0, -2.0, -3.0, -1.0, -1.0, -2.0, -2.0, -1.0, // C
	-1.0,  1.0,  0.0,  0.0, -3.0, 11.0,  2.0, -2.0,  0.0, -3.0, -2.0,  1.0,  0.0, -3.0, -1.0,  0.0, -1.0, -2.0, -1.0, -2.0, // Q
	-1.0,  0.0,  0.0,  2.0, -4.0,  2.0, 11.0, -2.0,  0.0, -3.0, -3.0,  1.0, -2.0, -3.0, -1.0,  0.0, -1.0, -3.0, -2.0, -2.0, // E
	 0.0, -2.0,  0.0, -1.0, -3.0, -2.0, -2.0, 11.0, -2.0, -4.0, -4.0, -2.0, -3.0, -3.0, -2.0,  0.0, -2.0, -2.0, -3.0, -3.0, // G
	-2.0,  0.0,  1.0, -1.0, -3.0,  0.0,  0.0, -2.0, 11.0, -3.0, -3.0, -1.0, -2.0, -1.0, -2.0, -1.0, -2.0, -2.0,  2.0, -3.0, // H
	-1.0, -3.0, -3.0, -3.0, -1.0, -3.0, -3.0, -4.0, -3.0, 11.0,  2.0, -3.0,  1.0,  0.0, -3.0, -2.0, -1.0, -3.0, -1.0,  3.0, // I
	-1.0, -2.0, -3.0, -4.0, -1.0, -2.0, -3.0, -4.0, -3.0,  2.0, 11.0, -2.0,  2.0,  0.0, -3.0, -2.0, -1.0, -2.0, -1.0,  1.0, // L
	-1.0,  2.0,  0.0, -1.0, -3.0,  1.0,  1.0, -2.0, -1.0, -3.0, -2.0, 11.0, -1.0, -3.0, -1.0,  0.0, -1.0, -3.0, -2.0, -2.0, // K
	-1.0, -1.0, -2.0, -3.0, -1.0,  0.0, -2.0, -3.0, -2.0,  1.0,  2.0, -1.0, 11.0,  0.0, -2.0, -1.0, -1.0, -1.0, -1.0,  1.0, // M
	-2.0, -3.0, -3.0, -3.0, -2.0, -3.0, -3.0, -3.0, -1.0,  0.0,  0.0, -3.0,  0.0, 11.0, -4.0, -2.0, -2.0,  1.0,  3.0, -1.0, // F
	-1.0, -2.0, -2.0, -1.0, -3.0, -1.0, -1.0, -2.0, -2.0, -3.0, -3.0, -1.0, -2.0, -4.0, 11.0, -1.0, -1.0, -4.0, -3.0, -2.0, // P
	 1.0, -1.0,  1.0,  0.0, -1.0,  0.0,  0.0,  0.0, -1.0, -2.0, -2.0,  0.0, -1.0, -2.0, -1.0, 11.0,  1.0, -3.0, -2.0, -2.0, // S
	 0.0, -1.0,  0.0, -1.0, -1.0, -1.0, -1.0, -2.0, -2.0, -1.0, -1.0, -1.0, -1.0, -2.0, -1.0,  1.0, 11.0, -2.0, -2.0,  0.0, // T
	-3.0, -3.0, -4.0, -4.0, -2.0, -2.0, -3.0, -2.0, -2.0, -3.0, -2.0, -3.0, -1.0,  1.0, -4.0, -3.0, -2.0, 11.0,  2.0, -3.0, // W
	-2.0, -2.0, -2.0, -3.0, -2.0, -1.0, -2.0, -3.0,  2.0, -1.0, -1.0, -2.0, -1.0,  3.0, -3.0, -2.0, -2.0,  2.0, 11.0, -1.0, // Y
	 0.0, -3.0, -3.0, -3.0, -1.0, -2.0, -2.0, -3.0, -3.0,  3.0,  1.0, -2.0,  1.0, -1.0, -2.0, -2.0,  0.0, -3.0, -1.0, 11.0  // V
	// A,    R,    N,    D,    C,    Q,    E,    G,    H,    I,    L,    K,    M,    F,    P,    S,    T,    W,    Y,    V
	];

	match ( *arg_m ).as_str() {
		"blosum45"    => _matrix = _blosum45,
		"blosum50"    => _matrix = _blosum50,
		"blosum62"    => _matrix = _blosum62,
		"blosum80"    => _matrix = _blosum80,
		"blosum90"    => _matrix = _blosum90,
		"pam30"       => _matrix = _pam30,
		"pam70"       => _matrix = _pam70,
		"pam250"      => _matrix = _pam250,
		"pet91mod"    => _matrix = _pet91mod,
		"blosum62mod" => _matrix = _blosum62mod,
		_             => _matrix = _blosum62,
	}

	_matrix
}