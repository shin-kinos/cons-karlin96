
pub fn define_matrix( arg_m : &String ) -> Vec<f64>
{
	/* 20 × 20 elements */
	let mut _matrix : Vec<f64> = vec![ 0.0; 400 ];

	let _blosum45 : Vec<f64> = vec![ // BLOSUM45
	   5., -2., -1., -2., -1., -1., -1.,  0., -2., -1., -1., -1., -1., -2., -1.,  1.,  0., -2., -2.,  0., // A
	  -2.,  7.,  0., -1., -3.,  1.,  0., -2.,  0., -3., -2.,  3., -1., -2., -2., -1., -1., -2., -1., -2., // R
	  -1.,  0.,  6.,  2., -2.,  0.,  0.,  0.,  1., -2., -3.,  0., -2., -2., -2.,  1.,  0., -4., -2., -3., // N
	  -2., -1.,  2.,  7., -3.,  0.,  2., -1.,  0., -4., -3.,  0., -3., -4., -1.,  0., -1., -4., -2., -3., // D
	  -1., -3., -2., -3., 12., -3., -3., -3., -3., -3., -2., -3., -2., -2., -4., -1., -1., -5., -3., -1., // C
	  -1.,  1.,  0.,  0., -3.,  6.,  2., -2.,  1., -2., -2.,  1.,  0., -4., -1.,  0., -1., -2., -1., -3., // Q
	  -1.,  0.,  0.,  2., -3.,  2.,  6., -2.,  0., -3., -2.,  1., -2., -3.,  0.,  0., -1., -3., -2., -3., // E
	   0., -2.,  0., -1., -3., -2., -2.,  7., -2., -4., -3., -2., -2., -3., -2.,  0., -2., -2., -3., -3., // G
	  -2.,  0.,  1.,  0., -3.,  1.,  0., -2., 10., -3., -2., -1.,  0., -2., -2., -1., -2., -3.,  2., -3., // H
	  -1., -3., -2., -4., -3., -2., -3., -4., -3.,  5.,  2., -3.,  2.,  0., -2., -2., -1., -2.,  0.,  3., // I
	  -1., -2., -3., -3., -2., -2., -2., -3., -2.,  2.,  5., -3.,  2.,  1., -3., -3., -1., -2.,  0.,  1., // L
	  -1.,  3.,  0.,  0., -3.,  1.,  1., -2., -1., -3., -3.,  5., -1., -3., -1., -1., -1., -2., -1., -2., // K
	  -1., -1., -2., -3., -2.,  0., -2., -2.,  0.,  2.,  2., -1.,  6.,  0., -2., -2., -1., -2.,  0.,  1., // M
	  -2., -2., -2., -4., -2., -4., -3., -3., -2.,  0.,  1., -3.,  0.,  8., -3., -2., -1.,  1.,  3.,  0., // F
	  -1., -2., -2., -1., -4., -1.,  0., -2., -2., -2., -3., -1., -2., -3.,  9., -1., -1., -3., -3., -3., // P
	   1., -1.,  1.,  0., -1.,  0.,  0.,  0., -1., -2., -3., -1., -2., -2., -1.,  4.,  2., -4., -2., -1., // S
	   0., -1.,  0., -1., -1., -1., -1., -2., -2., -1., -1., -1., -1., -1., -1.,  2.,  5., -3., -1.,  0., // T
	  -2., -2., -4., -4., -5., -2., -3., -2., -3., -2., -2., -2., -2.,  1., -3., -4., -3., 15.,  3., -3., // W
	  -2., -1., -2., -2., -3., -1., -2., -3.,  2.,  0.,  0., -1.,  0.,  3., -3., -2., -1.,  3.,  8., -1., // Y
	   0., -2., -3., -3., -1., -3., -3., -3., -3.,  3.,  1., -2.,  1.,  0., -3., -1.,  0., -3., -1.,  5.  // V
	// A    R    N    D    C    Q    E    G    H    I    L    K    M    F    P    S    T    W    Y    V
	];

	let _blosum50 : Vec<f64> = vec![ // BLOSUM50
	   5., -2., -1., -2., -1., -1., -1.,  0., -2., -1., -2., -1., -1., -3., -1.,  1.,  0., -3., -2.,  0., // A
	  -2.,  7., -1., -2., -4.,  1.,  0., -3.,  0., -4., -3.,  3., -2., -3., -3., -1., -1., -3., -1., -3., // R
	  -1., -1.,  7.,  2., -2.,  0.,  0.,  0.,  1., -3., -4.,  0., -2., -4., -2.,  1.,  0., -4., -2., -3., // N
	  -2., -2.,  2.,  8., -4.,  0.,  2., -1., -1., -4., -4., -1., -4., -5., -1.,  0., -1., -5., -3., -4., // D
	  -1., -4., -2., -4., 13., -3., -3., -3., -3., -2., -2., -3., -2., -2., -4., -1., -1., -5., -3., -1., // C
	  -1.,  1.,  0.,  0., -3.,  7.,  2., -2.,  1., -3., -2.,  2.,  0., -4., -1.,  0., -1., -1., -1., -3., // Q
	  -1.,  0.,  0.,  2., -3.,  2.,  6., -3.,  0., -4., -3.,  1., -2., -3., -1., -1., -1., -3., -2., -3., // E
	   0., -3.,  0., -1., -3., -2., -3.,  8., -2., -4., -4., -2., -3., -4., -2.,  0., -2., -3., -3., -4., // G
	  -2.,  0.,  1., -1., -3.,  1.,  0., -2., 10., -4., -3.,  0., -1., -1., -2., -1., -2., -3.,  2., -4., // H
	  -1., -4., -3., -4., -2., -3., -4., -4., -4.,  5.,  2., -3.,  2.,  0., -3., -3., -1., -3., -1.,  4., // I
	  -2., -3., -4., -4., -2., -2., -3., -4., -3.,  2.,  5., -3.,  3.,  1., -4., -3., -1., -2., -1.,  1., // L
	  -1.,  3.,  0., -1., -3.,  2.,  1., -2.,  0., -3., -3.,  6., -2., -4., -1.,  0., -1., -3., -2., -3., // K
	  -1., -2., -2., -4., -2.,  0., -2., -3., -1.,  2.,  3., -2.,  7.,  0., -3., -2., -1., -1.,  0.,  1., // M
	  -3., -3., -4., -5., -2., -4., -3., -4., -1.,  0.,  1., -4.,  0.,  8., -4., -3., -2.,  1.,  4., -1., // F
	  -1., -3., -2., -1., -4., -1., -1., -2., -2., -3., -4., -1., -3., -4., 10., -1., -1., -4., -3., -3., // P
	   1., -1.,  1.,  0., -1.,  0., -1.,  0., -1., -3., -3.,  0., -2., -3., -1.,  5.,  2., -4., -2., -2., // S
	   0., -1.,  0., -1., -1., -1., -1., -2., -2., -1., -1., -1., -1., -2., -1.,  2.,  5., -3., -2.,  0., // T
	  -3., -3., -4., -5., -5., -1., -3., -3., -3., -3., -2., -3., -1.,  1., -4., -4., -3., 15.,  2., -3., // W
	  -2., -1., -2., -3., -3., -1., -2., -3.,  2., -1., -1., -2.,  0.,  4., -3., -2., -2.,  2.,  8., -1., // Y
	   0., -3., -3., -4., -1., -3., -3., -4., -4.,  4.,  1., -3.,  1., -1., -3., -2.,  0., -3., -1.,  5.  // V
	// A    R    N    D    C    Q    E    G    H    I    L    K    M    F    P    S    T    W    Y    V
	];

	let _blosum62 : Vec<f64> = vec![ // BLOSUM62
	   4., -1., -2., -2.,  0., -1., -1.,  0., -2., -1., -1., -1., -1., -2., -1.,  1.,  0., -3., -2.,  0., // A
	  -1.,  5.,  0., -2., -3.,  1.,  0., -2.,  0., -3., -2.,  2., -1., -3., -2., -1., -1., -3., -2., -3., // R
	  -2.,  0.,  6.,  1., -3.,  0.,  0.,  0.,  1., -3., -3.,  0., -2., -3., -2.,  1.,  0., -4., -2., -3., // N
	  -2., -2.,  1.,  6., -3.,  0.,  2., -1., -1., -3., -4., -1., -3., -3., -1.,  0., -1., -4., -3., -3., // D
	   0., -3., -3., -3.,  9., -3., -4., -3., -3., -1., -1., -3., -1., -2., -3., -1., -1., -2., -2., -1., // C
	  -1.,  1.,  0.,  0., -3.,  5.,  2., -2.,  0., -3., -2.,  1.,  0., -3., -1.,  0., -1., -2., -1., -2., // Q
	  -1.,  0.,  0.,  2., -4.,  2.,  5., -2.,  0., -3., -3.,  1., -2., -3., -1.,  0., -1., -3., -2., -2., // E
	   0., -2.,  0., -1., -3., -2., -2.,  6., -2., -4., -4., -2., -3., -3., -2.,  0., -2., -2., -3., -3., // G
	  -2.,  0.,  1., -1., -3.,  0.,  0., -2.,  8., -3., -3., -1., -2., -1., -2., -1., -2., -2.,  2., -3., // H
	  -1., -3., -3., -3., -1., -3., -3., -4., -3.,  4.,  2., -3.,  1.,  0., -3., -2., -1., -3., -1.,  3., // I
	  -1., -2., -3., -4., -1., -2., -3., -4., -3.,  2.,  4., -2.,  2.,  0., -3., -2., -1., -2., -1.,  1., // L
	  -1.,  2.,  0., -1., -3.,  1.,  1., -2., -1., -3., -2.,  5., -1., -3., -1.,  0., -1., -3., -2., -2., // K
	  -1., -1., -2., -3., -1.,  0., -2., -3., -2.,  1.,  2., -1.,  5.,  0., -2., -1., -1., -1., -1.,  1., // M
	  -2., -3., -3., -3., -2., -3., -3., -3., -1.,  0.,  0., -3.,  0.,  6., -4., -2., -2.,  1.,  3., -1., // F
	  -1., -2., -2., -1., -3., -1., -1., -2., -2., -3., -3., -1., -2., -4.,  7., -1., -1., -4., -3., -2., // P
	   1., -1.,  1.,  0., -1.,  0.,  0.,  0., -1., -2., -2.,  0., -1., -2., -1.,  4.,  1., -3., -2., -2., // S
	   0., -1.,  0., -1., -1., -1., -1., -2., -2., -1., -1., -1., -1., -2., -1.,  1.,  5., -2., -2.,  0., // T
	  -3., -3., -4., -4., -2., -2., -3., -2., -2., -3., -2., -3., -1.,  1., -4., -3., -2., 11.,  2., -3., // W
	  -2., -2., -2., -3., -2., -1., -2., -3.,  2., -1., -1., -2., -1.,  3., -3., -2., -2.,  2.,  7., -1., // Y
	   0., -3., -3., -3., -1., -2., -2., -3., -3.,  3.,  1., -2.,  1., -1., -2., -2.,  0., -3., -1.,  4.  // V
	// A    R    N    D    C    Q    E    G    H    I    L    K    M    F    P    S    T    W    Y    V
	];

	let _blosum80 : Vec<f64> = vec![ // BLOSUM80
	   5., -2., -2., -2., -1., -1., -1.,  0., -2., -2., -2., -1., -1., -3., -1.,  1.,  0., -3., -2.,  0., // A
	  -2.,  6., -1., -2., -4.,  1., -1., -3.,  0., -3., -3.,  2., -2., -4., -2., -1., -1., -4., -3., -3., // R
	  -2., -1.,  6.,  1., -3.,  0., -1., -1.,  0., -4., -4.,  0., -3., -4., -3.,  0.,  0., -4., -3., -4., // N
	  -2., -2.,  1.,  6., -4., -1.,  1., -2., -2., -4., -5., -1., -4., -4., -2., -1., -1., -6., -4., -4., // D
	  -1., -4., -3., -4.,  9., -4., -5., -4., -4., -2., -2., -4., -2., -3., -4., -2., -1., -3., -3., -1., // C
	  -1.,  1.,  0., -1., -4.,  6.,  2., -2.,  1., -3., -3.,  1.,  0., -4., -2.,  0., -1., -3., -2., -3., // Q
	  -1., -1., -1.,  1., -5.,  2.,  6., -3.,  0., -4., -4.,  1., -2., -4., -2.,  0., -1., -4., -3., -3., // E
	   0., -3., -1., -2., -4., -2., -3.,  6., -3., -5., -4., -2., -4., -4., -3., -1., -2., -4., -4., -4., // G
	  -2.,  0.,  0., -2., -4.,  1.,  0., -3.,  8., -4., -3., -1., -2., -2., -3., -1., -2., -3.,  2., -4., // H
	  -2., -3., -4., -4., -2., -3., -4., -5., -4.,  5.,  1., -3.,  1., -1., -4., -3., -1., -3., -2.,  3., // I
	  -2., -3., -4., -5., -2., -3., -4., -4., -3.,  1.,  4., -3.,  2.,  0., -3., -3., -2., -2., -2.,  1., // L
	  -1.,  2.,  0., -1., -4.,  1.,  1., -2., -1., -3., -3.,  5., -2., -4., -1., -1., -1., -4., -3., -3., // K
	  -1., -2., -3., -4., -2.,  0., -2., -4., -2.,  1.,  2., -2.,  6.,  0., -3., -2., -1., -2., -2.,  1., // M
	  -3., -4., -4., -4., -3., -4., -4., -4., -2., -1.,  0., -4.,  0.,  6., -4., -3., -2.,  0.,  3., -1., // F
	  -1., -2., -3., -2., -4., -2., -2., -3., -3., -4., -3., -1., -3., -4.,  8., -1., -2., -5., -4., -3., // P
	   1., -1.,  0., -1., -2.,  0.,  0., -1., -1., -3., -3., -1., -2., -3., -1.,  5.,  1., -4., -2., -2., // S
	   0., -1.,  0., -1., -1., -1., -1., -2., -2., -1., -2., -1., -1., -2., -2.,  1.,  5., -4., -2.,  0., // T
	  -3., -4., -4., -6., -3., -3., -4., -4., -3., -3., -2., -4., -2.,  0., -5., -4., -4., 11.,  2., -3., // W
	  -2., -3., -3., -4., -3., -2., -3., -4.,  2., -2., -2., -3., -2.,  3., -4., -2., -2.,  2.,  7., -2., // Y
	   0., -3., -4., -4., -1., -3., -3., -4., -4.,  3.,  1., -3.,  1., -1., -3., -2.,  0., -3., -2.,  4.  // V
	// A    R    N    D    C    Q    E    G    H    I    L    K    M    F    P    S    T    W    Y    V
	];

	let _blosum90 : Vec<f64> = vec![ // BLOSUM90
	   5., -2., -2., -3., -1., -1., -1.,  0., -2., -2., -2., -1., -2., -3., -1.,  1.,  0., -4., -3., -1., // A
	  -2.,  6., -1., -3., -5.,  1., -1., -3.,  0., -4., -3.,  2., -2., -4., -3., -1., -2., -4., -3., -3., // R
	  -2., -1.,  7.,  1., -4.,  0., -1., -1.,  0., -4., -4.,  0., -3., -4., -3.,  0.,  0., -5., -3., -4., // N
	  -3., -3.,  1.,  7., -5., -1.,  1., -2., -2., -5., -5., -1., -4., -5., -3., -1., -2., -6., -4., -5., // D
	  -1., -5., -4., -5.,  9., -4., -6., -4., -5., -2., -2., -4., -2., -3., -4., -2., -2., -4., -4., -2., // C
	  -1.,  1.,  0., -1., -4.,  7.,  2., -3.,  1., -4., -3.,  1.,  0., -4., -2., -1., -1., -3., -3., -3., // Q
	  -1., -1., -1.,  1., -6.,  2.,  6., -3., -1., -4., -4.,  0., -3., -5., -2., -1., -1., -5., -4., -3., // E
	   0., -3., -1., -2., -4., -3., -3.,  6., -3., -5., -5., -2., -4., -5., -3., -1., -3., -4., -5., -5., // G
	  -2.,  0.,  0., -2., -5.,  1., -1., -3.,  8., -4., -4., -1., -3., -2., -3., -2., -2., -3.,  1., -4., // H
	  -2., -4., -4., -5., -2., -4., -4., -5., -4.,  5.,  1., -4.,  1., -1., -4., -3., -1., -4., -2.,  3., // I
	  -2., -3., -4., -5., -2., -3., -4., -5., -4.,  1.,  5., -3.,  2.,  0., -4., -3., -2., -3., -2.,  0., // L
	  -1.,  2.,  0., -1., -4.,  1.,  0., -2., -1., -4., -3.,  6., -2., -4., -2., -1., -1., -5., -3., -3., // K
	  -2., -2., -3., -4., -2.,  0., -3., -4., -3.,  1.,  2., -2.,  7., -1., -3., -2., -1., -2., -2.,  0., // M
	  -3., -4., -4., -5., -3., -4., -5., -5., -2., -1.,  0., -4., -1.,  7., -4., -3., -3.,  0.,  3., -2., // F
	  -1., -3., -3., -3., -4., -2., -2., -3., -3., -4., -4., -2., -3., -4.,  8., -2., -2., -5., -4., -3., // P
	   1., -1.,  0., -1., -2., -1., -1., -1., -2., -3., -3., -1., -2., -3., -2.,  5.,  1., -4., -3., -2., // S
	   0., -2.,  0., -2., -2., -1., -1., -3., -2., -1., -2., -1., -1., -3., -2.,  1.,  6., -4., -2., -1., // T
	  -4., -4., -5., -6., -4., -3., -5., -4., -3., -4., -3., -5., -2.,  0., -5., -4., -4., 11.,  2., -3., // W
	  -3., -3., -3., -4., -4., -3., -4., -5.,  1., -2., -2., -3., -2.,  3., -4., -3., -2.,  2.,  8., -3., // Y
	  -1., -3., -4., -5., -2., -3., -3., -5., -4.,  3.,  0., -3.,  0., -2., -3., -2., -1., -3., -3.,  5.  // V
	// A    R    N    D    C    Q    E    G    H    I    L    K    M    F    P    S    T    W    Y    V
	];

	let _pam30 : Vec<f64> = vec![ // PAM30
	   6., -7., -4., -3., -6., -4., -2., -2., -7., -5., -6., -7., -5., -8., -2.,  0., -1.,-13., -8., -2., // A
	  -7.,  8., -6.,-10., -8., -2., -9., -9., -2., -5., -8.,  0., -4., -9., -4., -3., -6., -2.,-10., -8., // R
	  -4., -6.,  8.,  2.,-11., -3., -2., -3.,  0., -5., -7., -1., -9., -9., -6.,  0., -2., -8., -4., -8., // N
	  -3.,-10.,  2.,  8.,-14., -2.,  2., -3., -4., -7.,-12., -4.,-11.,-15., -8., -4., -5.,-15.,-11., -8., // D
	  -6., -8.,-11.,-14., 10.,-14.,-14., -9., -7., -6.,-15.,-14.,-13.,-13., -8., -3., -8.,-15., -4., -6., // C
	  -4., -2., -3., -2.,-14.,  8.,  1., -7.,  1., -8., -5., -3., -4.,-13., -3., -5., -5.,-13.,-12., -7., // Q
	  -2., -9., -2.,  2.,-14.,  1.,  8., -4., -5., -5., -9., -4., -7.,-14., -5., -4., -6.,-17., -8., -6., // E
	  -2., -9., -3., -3., -9., -7., -4.,  6., -9.,-11.,-10., -7., -8., -9., -6., -2., -6.,-15.,-14., -5., // G
	  -7., -2.,  0., -4., -7.,  1., -5., -9.,  9., -9., -6., -6.,-10., -6., -4., -6., -7., -7., -3., -6., // H
	  -5., -5., -5., -7., -6., -8., -5.,-11., -9.,  8., -1., -6., -1., -2., -8., -7., -2.,-14., -6.,  2., // I
	  -6., -8., -7.,-12.,-15., -5., -9.,-10., -6., -1.,  7., -8.,  1., -3., -7., -8., -7., -6., -7., -2., // L
	  -7.,  0., -1., -4.,-14., -3., -4., -7., -6., -6., -8.,  7., -2.,-14., -6., -4., -3.,-12., -9., -9., // K
	  -5., -4., -9.,-11.,-13., -4., -7., -8.,-10., -1.,  1., -2., 11., -4., -8., -5., -4.,-13.,-11., -1., // M
	  -8., -9., -9.,-15.,-13.,-13.,-14., -9., -6., -2., -3.,-14., -4.,  9.,-10., -6., -9., -4.,  2., -8., // F
	  -2., -4., -6., -8., -8., -3., -5., -6., -4., -8., -7., -6., -8.,-10.,  8., -2., -4.,-14.,-13., -6., // P
	   0., -3.,  0., -4., -3., -5., -4., -2., -6., -7., -8., -4., -5., -6., -2.,  6.,  0., -5., -7., -6., // S
	  -1., -6., -2., -5., -8., -5., -6., -6., -7., -2., -7., -3., -4., -9., -4.,  0.,  7.,-13., -6., -3., // T
	 -13., -2., -8.,-15.,-15.,-13.,-17.,-15., -7.,-14., -6.,-12.,-13., -4.,-14., -5.,-13., 13., -5.,-15., // W
	  -8.,-10., -4.,-11., -4.,-12., -8.,-14., -3., -6., -7., -9.,-11.,  2.,-13., -7., -6., -5., 10., -7., // Y
	  -2., -8., -8., -8., -6., -7., -6., -5., -6.,  2., -2., -9., -1., -8., -6., -6., -3.,-15., -7.,  7.  // V
	// A    R    N    D    C    Q    E    G    H    I    L    K    M    F    P    S    T    W    Y    V
	];

	let _pam70 : Vec<f64> = vec![ // PAM70
	   5., -4., -2., -1., -4., -2., -1.,  0., -4., -2., -4., -4., -3., -6.,  0.,  1.,  1., -9., -5., -1., // A
	  -4.,  8., -3., -6., -5.,  0., -5., -6.,  0., -3., -6.,  2., -2., -7., -2., -1., -4.,  0., -7., -5., // R
	  -2., -3.,  6.,  3., -7., -1.,  0., -1.,  1., -3., -5.,  0., -5., -6., -3.,  1.,  0., -6., -3., -5., // N
	  -1., -6.,  3.,  6., -9.,  0.,  3., -1., -1., -5., -8., -2., -7.,-10., -4., -1., -2.,-10., -7., -5., // D
	  -4., -5., -7., -9.,  9., -9., -9., -6., -5., -4.,-10., -9., -9., -8., -5., -1., -5.,-11., -2., -4., // C
	  -2.,  0., -1.,  0., -9.,  7.,  2., -4.,  2., -5., -3., -1., -2., -9., -1., -3., -3., -8., -8., -4., // Q
	  -1., -5.,  0.,  3., -9.,  2.,  6., -2., -2., -4., -6., -2., -4., -9., -3., -2., -3.,-11., -6., -4., // E
	   0., -6., -1., -1., -6., -4., -2.,  6., -6., -6., -7., -5., -6., -7., -3.,  0., -3.,-10., -9., -3., // G
	  -4.,  0.,  1., -1., -5.,  2., -2., -6.,  8., -6., -4., -3., -6., -4., -2., -3., -4., -5., -1., -4., // H
	  -2., -3., -3., -5., -4., -5., -4., -6., -6.,  7.,  1., -4.,  1.,  0., -5., -4., -1., -9., -4.,  3., // I
	  -4., -6., -5., -8.,-10., -3., -6., -7., -4.,  1.,  6., -5.,  2., -1., -5., -6., -4., -4., -4.,  0., // L
	  -4.,  2.,  0., -2., -9., -1., -2., -5., -3., -4., -5.,  6.,  0., -9., -4., -2., -1., -7., -7., -6., // K
	  -3., -2., -5., -7., -9., -2., -4., -6., -6.,  1.,  2.,  0., 10., -2., -5., -3., -2., -8., -7.,  0., // M
	  -6., -7., -6.,-10., -8., -9., -9., -7., -4.,  0., -1., -9., -2.,  8., -7., -4., -6., -2.,  4., -5., // F
	   0., -2., -3., -4., -5., -1., -3., -3., -2., -5., -5., -4., -5., -7.,  7.,  0., -2., -9., -9., -3., // P
	   1., -1.,  1., -1., -1., -3., -2.,  0., -3., -4., -6., -2., -3., -4.,  0.,  5.,  2., -3., -5., -3., // S
	   1., -4.,  0., -2., -5., -3., -3., -3., -4., -1., -4., -1., -2., -6., -2.,  2.,  6., -8., -4., -1., // T
	  -9.,  0., -6.,-10.,-11., -8.,-11.,-10., -5., -9., -4., -7., -8., -2., -9., -3., -8., 13., -3.,-10., // W
	  -5., -7., -3., -7., -2., -8., -6., -9., -1., -4., -4., -7., -7.,  4., -9., -5., -4., -3.,  9., -5., // Y
	  -1., -5., -5., -5., -4., -4., -4., -3., -4.,  3.,  0., -6.,  0., -5., -3., -3., -1.,-10., -5.,  6.  // V
	// A    R    N    D    C    Q    E    G    H    I    L    K    M    F    P    S    T    W    Y    V
	];

	let _pam250 : Vec<f64> = vec![ // PAM250
	   2., -2.,  0.,  0., -2.,  0.,  0.,  1., -1., -1., -2., -1., -1., -3.,  1.,  1.,  1., -6., -3.,  0., // A
	  -2.,  6.,  0., -1., -4.,  1., -1., -3.,  2., -2., -3.,  3.,  0., -4.,  0.,  0., -1.,  2., -4., -2., // R
	   0.,  0.,  2.,  2., -4.,  1.,  1.,  0.,  2., -2., -3.,  1., -2., -3.,  0.,  1.,  0., -4., -2., -2., // N
	   0., -1.,  2.,  4., -5.,  2.,  3.,  1.,  1., -2., -4.,  0., -3., -6., -1.,  0.,  0., -7., -4., -2., // D
	  -2., -4., -4., -5., 12., -5., -5., -3., -3., -2., -6., -5., -5., -4., -3.,  0., -2., -8.,  0., -2., // C
	   0.,  1.,  1.,  2., -5.,  4.,  2., -1.,  3., -2., -2.,  1., -1., -5.,  0., -1., -1., -5., -4., -2., // Q
	   0., -1.,  1.,  3., -5.,  2.,  4.,  0.,  1., -2., -3.,  0., -2., -5., -1.,  0.,  0., -7., -4., -2., // E
	   1., -3.,  0.,  1., -3., -1.,  0.,  5., -2., -3., -4., -2., -3., -5.,  0.,  1.,  0., -7., -5., -1., // G
	  -1.,  2.,  2.,  1., -3.,  3.,  1., -2.,  6., -2., -2.,  0., -2., -2.,  0., -1., -1., -3.,  0., -2., // H
	  -1., -2., -2., -2., -2., -2., -2., -3., -2.,  5.,  2., -2.,  2.,  1., -2., -1.,  0., -5., -1.,  4., // I
	  -2., -3., -3., -4., -6., -2., -3., -4., -2.,  2.,  6., -3.,  4.,  2., -3., -3., -2., -2., -1.,  2., // L
	  -1.,  3.,  1.,  0., -5.,  1.,  0., -2.,  0., -2., -3.,  5.,  0., -5., -1.,  0.,  0., -3., -4., -2., // K
	  -1.,  0., -2., -3., -5., -1., -2., -3., -2.,  2.,  4.,  0.,  6.,  0., -2., -2., -1., -4., -2.,  2., // M
	  -3., -4., -3., -6., -4., -5., -5., -5., -2.,  1.,  2., -5.,  0.,  9., -5., -3., -3.,  0.,  7., -1., // F
	   1.,  0.,  0., -1., -3.,  0., -1.,  0.,  0., -2., -3., -1., -2., -5.,  6.,  1.,  0., -6., -5., -1., // P
	   1.,  0.,  1.,  0.,  0., -1.,  0.,  1., -1., -1., -3.,  0., -2., -3.,  1.,  2.,  1., -2., -3., -1., // S
	   1., -1.,  0.,  0., -2., -1.,  0.,  0., -1.,  0., -2.,  0., -1., -3.,  0.,  1.,  3., -5., -3.,  0., // T
	  -6.,  2., -4., -7., -8., -5., -7., -7., -3., -5., -2., -3., -4.,  0., -6., -2., -5., 17.,  0., -6., // W
	  -3., -4., -2., -4.,  0., -4., -4., -5.,  0., -1., -1., -4., -2.,  7., -5., -3., -3.,  0., 10., -2., // Y
	   0., -2., -2., -2., -2., -2., -2., -1., -2.,  4.,  2., -2.,  2., -1., -1., -1.,  0., -6., -2.,  4.  // V
	// A    R    N    D    C    Q    E    G    H    I    L    K    M    F    P    S    T    W    Y    V
	];

	let _pet91mod : Vec<f64> = vec![ // Modified version of PET91
	  15., -1.,  0., -1., -1., -1., -1.,  1., -2.,  0., -1., -1., -1., -3.,  1.,  1.,  2., -4., -3.,  1., // A
	  -1., 15.,  0., -1., -1.,  2.,  0.,  0.,  2., -3., -3.,  4., -2., -4., -1., -1., -1.,  0., -2., -3., // R
	   0.,  0., 15.,  2., -1.,  0.,  1.,  0.,  1., -2., -3.,  1., -2., -3., -1.,  1.,  1., -4., -1., -2., // N
	  -1., -1.,  2., 15., -3.,  0.,  4.,  1.,  0., -3., -4.,  0., -3., -5., -2.,  0., -1., -5., -2., -3., // D
	  -1., -1., -1., -3., 15., -3., -4., -1.,  0., -2., -3., -3., -2.,  0., -2.,  1., -1.,  1.,  2., -2., // C
	  -1.,  2.,  0.,  0., -3., 15.,  2., -1.,  3., -3., -2.,  2., -2., -4.,  0., -1., -1., -3., -1., -3., // Q
	  -1.,  0.,  1.,  4., -4.,  2., 15.,  1.,  0., -3., -4.,  1., -3., -5., -2., -1., -1., -5., -4., -2., // E
	   1.,  0.,  0.,  1., -1., -1.,  1., 15., -2., -3., -4., -1., -3., -5., -1.,  1.,  0., -2., -4., -2., // G
	  -2.,  2.,  1.,  0.,  0.,  3.,  0., -2., 15., -3., -2.,  1., -2.,  0.,  0., -1., -1., -3.,  4., -3., // H
	   0., -3., -2., -3., -2., -3., -3., -3., -3., 15.,  2., -3.,  3.,  0., -2., -1.,  1., -4., -2.,  4., // I
	  -1., -3., -3., -4., -3., -2., -4., -4., -2.,  2., 15., -3.,  3.,  2.,  0., -2., -1., -2., -1.,  2., // L
	  -1.,  4.,  1.,  0., -3.,  2.,  1., -1.,  1., -3., -3., 15., -2., -5., -2., -1., -1., -3., -3., -3., // K
	  -1., -2., -2., -3., -2., -2., -3., -3., -2.,  3.,  3., -2., 15.,  0., -2., -1.,  0., -3., -3.,  2., // M
	  -3., -4., -3., -5.,  0., -4., -5., -5.,  0.,  0.,  2., -5.,  0., 15., -2., -2., -2., -1.,  5.,  0., // F
	   1., -1., -1., -2., -2.,  0., -2., -1.,  0., -2.,  0., -2., -2., -2., 15.,  1.,  1., -5., -3., -1., // P
	   1., -1.,  1.,  0.,  1., -1., -1.,  1., -1., -1., -2., -1., -1., -2.,  1., 15.,  1., -3., -1., -1., // S
	   2., -1.,  1., -1., -1., -1., -1.,  0., -1.,  1., -1., -1.,  0., -2.,  1.,  1., 15., -4., -3.,  0., // T
	  -4.,  0., -4., -5.,  1., -3., -5., -2., -3., -4., -2., -3., -3., -1., -5., -3., -4., 15.,  0., -4., // W
	  -3., -2., -1., -2.,  2., -1., -4., -4.,  4., -2., -1., -3., -3.,  5., -3., -1., -3.,  0., 15., -3., // Y
	   1., -3., -2., -3., -2., -3., -2., -2., -3.,  4.,  2., -3.,  2.,  0., -1., -1.,  0., -4., -3., 15.  // V
	// A    R    N    D    C    Q    E    G    H    I    L    K    M    F    P    S    T    W    Y    V
	];

	let _blosum62mod : Vec<f64> = vec![ // Modified version of BLOSUM62
	  11., -1., -2., -2.,  0., -1., -1.,  0., -2., -1., -1., -1., -1., -2., -1.,  1.,  0., -3., -2.,  0., // A
	  -1., 11.,  0., -2., -3.,  1.,  0., -2.,  0., -3., -2.,  2., -1., -3., -2., -1., -1., -3., -2., -3., // R
	  -2.,  0., 11.,  1., -3.,  0.,  0.,  0.,  1., -3., -3.,  0., -2., -3., -2.,  1.,  0., -4., -2., -3., // N
	  -2., -2.,  1., 11., -3.,  0.,  2., -1., -1., -3., -4., -1., -3., -3., -1.,  0., -1., -4., -3., -3., // D
	   0., -3., -3., -3., 11., -3., -4., -3., -3., -1., -1., -3., -1., -2., -3., -1., -1., -2., -2., -1., // C
	  -1.,  1.,  0.,  0., -3., 11.,  2., -2.,  0., -3., -2.,  1.,  0., -3., -1.,  0., -1., -2., -1., -2., // Q
	  -1.,  0.,  0.,  2., -4.,  2., 11., -2.,  0., -3., -3.,  1., -2., -3., -1.,  0., -1., -3., -2., -2., // E
	   0., -2.,  0., -1., -3., -2., -2., 11., -2., -4., -4., -2., -3., -3., -2.,  0., -2., -2., -3., -3., // G
	  -2.,  0.,  1., -1., -3.,  0.,  0., -2., 11., -3., -3., -1., -2., -1., -2., -1., -2., -2.,  2., -3., // H
	  -1., -3., -3., -3., -1., -3., -3., -4., -3., 11.,  2., -3.,  1.,  0., -3., -2., -1., -3., -1.,  3., // I
	  -1., -2., -3., -4., -1., -2., -3., -4., -3.,  2., 11., -2.,  2.,  0., -3., -2., -1., -2., -1.,  1., // L
	  -1.,  2.,  0., -1., -3.,  1.,  1., -2., -1., -3., -2., 11., -1., -3., -1.,  0., -1., -3., -2., -2., // K
	  -1., -1., -2., -3., -1.,  0., -2., -3., -2.,  1.,  2., -1., 11.,  0., -2., -1., -1., -1., -1.,  1., // M
	  -2., -3., -3., -3., -2., -3., -3., -3., -1.,  0.,  0., -3.,  0., 11., -4., -2., -2.,  1.,  3., -1., // F
	  -1., -2., -2., -1., -3., -1., -1., -2., -2., -3., -3., -1., -2., -4., 11., -1., -1., -4., -3., -2., // P
	   1., -1.,  1.,  0., -1.,  0.,  0.,  0., -1., -2., -2.,  0., -1., -2., -1., 11.,  1., -3., -2., -2., // S
	   0., -1.,  0., -1., -1., -1., -1., -2., -2., -1., -1., -1., -1., -2., -1.,  1., 11., -2., -2.,  0., // T
	  -3., -3., -4., -4., -2., -2., -3., -2., -2., -3., -2., -3., -1.,  1., -4., -3., -2., 11.,  2., -3., // W
	  -2., -2., -2., -3., -2., -1., -2., -3.,  2., -1., -1., -2., -1.,  3., -3., -2., -2.,  2., 11., -1., // Y
	   0., -3., -3., -3., -1., -2., -2., -3., -3.,  3.,  1., -2.,  1., -1., -2., -2.,  0., -3., -1., 11.  // V
	// A    R    N    D    C    Q    E    G    H    I    L    K    M    F    P    S    T    W    Y    V
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
