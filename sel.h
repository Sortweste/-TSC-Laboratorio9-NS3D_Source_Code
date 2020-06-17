node selectNode(int i, element e,mesh &m){
	node n;
	switch(i){
		case 1: n = m.getNode(e.getNode1()-1); break;
		case 2: n = m.getNode(e.getNode2()-1); break;
		case 3: n = m.getNode(e.getNode3()-1); break;
	}
	return n;
}

float selectCoord(int c, node n){
	float v;
	switch(c){
		case EQUIS: v = n.getX(); break;
		case YE: v = n.getY(); break;
	}
	return v;
}

float calcularTenedor(element e, int coord, int i, int j,mesh &m){
	node n1=selectNode(i,e,m),n2=selectNode(j,e,m);

	return selectCoord(coord,n1) - selectCoord(coord,n2);
}

float calculateLocalD(int i,mesh m){
    Matrix matriz;
    Vector row1, row2;

    element e = m.getElement(i);
    node n1 = m.getNode(e.getNode1()-1);
    node n2 = m.getNode(e.getNode2()-1);
    node n3 = m.getNode(e.getNode3()-1);

	row1.push_back(calcularTenedor(e,EQUIS,2,1,m)); row1.push_back(calcularTenedor(e,YE,2,1,m));
	row2.push_back(calcularTenedor(e,EQUIS,3,1,m)); row2.push_back(calcularTenedor(e,YE,3,1,m));

	matriz.push_back(row1); matriz.push_back(row2);

    return determinant(matriz);
}

float calculateMagnitude(float v1, float v2){
    return sqrt(pow(v1,2)+pow(v2,2));
}

float calculateLocalArea(int i,mesh m){
    float A,s,a,b,c;
    element e = m.getElement(i);
    node n1 = m.getNode(e.getNode1()-1);
    node n2 = m.getNode(e.getNode2()-1);
    node n3 = m.getNode(e.getNode3()-1);

    a = calculateMagnitude(n2.getX()-n1.getX(),n2.getY()-n1.getY());
    b = calculateMagnitude(n3.getX()-n2.getX(),n3.getY()-n2.getY());
    c = calculateMagnitude(n3.getX()-n1.getX(),n3.getY()-n1.getY());
    s = (a+b+c)/2;

    A = sqrt(s*(s-a)*(s-b)*(s-c));
    return A;
}

void calculateLocalA(int i,Matrix &A,mesh m){
    zeroes(A,2);
    element e = m.getElement(i);
    node n1 = m.getNode(e.getNode1()-1);
    node n2 = m.getNode(e.getNode2()-1);
    node n3 = m.getNode(e.getNode3()-1);

    A.at(0).at(0) = calcularTenedor(e,YE,3,1,m);  A.at(0).at(1) = calcularTenedor(e,YE,1,2,m);
    A.at(1).at(0) = calcularTenedor(e,EQUIS,1,3,m);  A.at(1).at(1) = calcularTenedor(e,EQUIS,2,1,m);
}

//Matriz Beta
void calculateBetaMatrix(Matrix &B){
    zeroes(B,2,6);
    B.at(0).at(0) = -1; B.at(0).at(1) = 1; B.at(0).at(2) = 0; B.at(0).at(3) = -1; B.at(0).at(4) = 1; B.at(0).at(5) = 0;
    B.at(1).at(0) = -1; B.at(1).at(1) = 0; B.at(1).at(2) = 1; B.at(1).at(3) = -1; B.at(1).at(4) = 0; B.at(1).at(5) = 1;
}

void calculateBPrima(Matrix &C){
    zeroes(C,2,3);
    C.at(0).at(0) = -1; C.at(0).at(1) = 1; C.at(0).at(2) = 0;
    C.at(1).at(0) = -1; C.at(1).at(1) = 0; C.at(1).at(2) = 1;
}

void ubicarSubMatriz(Matrix &K,int fi,int ff,int ci,int cf,Matrix M){
    int n = 0, m= 0;
    for(int i=fi;i<=ff;i++){
        for(int j=ci;j<=cf;j++){
            K.at(i).at(j) = M.at(n).at(m);
            m++;
        }
        n++; m = 0;
    }
}

void calculateGammaMatrix(Matrix& m){
	zeroes(m,6,2);

	m.at(0).at(0) = 1;   m.at(0).at(1) = 0;
	m.at(1).at(0) = 1; m.at(1).at(1) = 0;
	m.at(2).at(0) = 1; m.at(2).at(1) = 0;
	m.at(3).at(0) = 0;   m.at(3).at(1) = 1;
	m.at(4).at(0) = 0;   m.at(4).at(1) = 1;
	m.at(5).at(0) = 0;   m.at(5).at(1) = 1;
}

float calculateLocalJ(int i,mesh m){
    Matrix matriz;
    Vector row1, row2;

    element e = m.getElement(i);
    node n1 = m.getNode(e.getNode1()-1);
    node n2 = m.getNode(e.getNode2()-1);
    node n3 = m.getNode(e.getNode3()-1);

	row1.push_back(calcularTenedor(e,EQUIS,2,1,m)); row1.push_back(calcularTenedor(e,EQUIS,3,1,m));
	row2.push_back(calcularTenedor(e,YE,2,1,m)); row2.push_back(calcularTenedor(e,YE,3,1,m));

	matriz.push_back(row1); matriz.push_back(row2);

    return determinant(matriz);
}

Matrix createLocalK(int e,mesh &m){
    //Preparaci�n de ingredientes
    float u_bar,nu,rho,Ae,J,D;
    
    //Componentes de K
    // [ A+K  G ]
    // [  D   0 ]
    Matrix matrixA,matrixK,matrixG,matrixD;
    Matrix K,g_matrix,g_matrix_t,Alpha,Beta,Alphat,Betat,BPrima,BPrimat;

    //Preparando matrixA (En clase conocida simplemente como A)
    u_bar = m.getParameter(ADJECTIVE_VELOCITY);
    J = calculateLocalJ(e,m);
    D = calculateLocalD(e,m);

    if(D == 0){
        cout << "\n!---CATASTROPHIC FAILURE---!\n";
        exit(EXIT_FAILURE);
    }
    
    calculateGammaMatrix(g_matrix);
    calculateLocalA(e,Alpha,m);
    calculateBetaMatrix(Beta);
    productRealMatrix(u_bar*J/(6*D),productMatrixMatrix(g_matrix,productMatrixMatrix(Alpha,Beta,2,2,6),6,2,6),matrixA);

    //Preparando matrixK (En clase conocida simplemente como K)
    nu = m.getParameter(DYNAMIC_VISCOSITY);
    Ae = calculateLocalArea(e,m);
    transpose(Alpha,Alphat);
    transpose(Beta,Betat);
    productRealMatrix(nu*Ae/(D*D),productMatrixMatrix(Betat,productMatrixMatrix(Alphat,productMatrixMatrix(Alpha,Beta,2,2,6),2,2,6),6,2,6),matrixK);

    //Preparando matrixG (En clase conocida simplemente como G)
    rho = m.getParameter(DENSITY);
    calculateBPrima(BPrima);
    productRealMatrix(J/(6*rho*D),productMatrixMatrix(g_matrix,productMatrixMatrix(Alpha,BPrima,2,2,3),6,2,3),matrixG);

    //Preparando matrixD (En clase conocida simplemente como D)
    transpose(BPrima,BPrimat);
    transpose(g_matrix,g_matrix_t);
    productRealMatrix(J/(6*D),productMatrixMatrix(BPrimat,productMatrixMatrix(Alphat,g_matrix_t,2,2,6),3,2,6),matrixD);

    //Colocando submatrices en K
    zeroes(K,9);
    ubicarSubMatriz(K,0,5,0,5,sumMatrix(matrixA,matrixK,6,6));
    ubicarSubMatriz(K,0,5,6,8,matrixG);
    ubicarSubMatriz(K,6,8,0,5,matrixD );

    return K;
}

Vector createLocalb(int e,mesh &m){
    Vector b0,b,f;
    Matrix g_matrix;

    float f_x = m.getParameter(EXTERNAL_FORCE_X);
    float f_y = m.getParameter(EXTERNAL_FORCE_Y);
    float J = calculateLocalJ(e,m);
    calculateGammaMatrix(g_matrix);
    zeroes(f,2);
    f.at(0) = f_x;
    f.at(1) = f_y;

    zeroes(b0,6);
    productMatrixVector(g_matrix,f,b0);
    productRealVector(J/6,b0,b);
    b.push_back(0); b.push_back(0); b.push_back(0);

    return b;
}

void crearSistemasLocales(mesh &m,vector<Matrix> &localKs,vector<Vector> &localbs){
    for(int i=0;i<m.getSize(ELEMENTS);i++){
        localKs.push_back(createLocalK(i,m));
        localbs.push_back(createLocalb(i,m));
    }
}

void assemblyK(element e,Matrix localK,Matrix &K,int nnodes){
    int index1 = e.getNode1() - 1;
    int index2 = e.getNode2() - 1;
    int index3 = e.getNode3() - 1;
    int index4 = index1+nnodes, index5 = index2+nnodes, index6 = index3+nnodes;
    int index7 = index1+2*nnodes, index8 = index2+2*nnodes, index9 = index3+2*nnodes;

    K.at(index1).at(index1) += localK.at(0).at(0);
    K.at(index1).at(index2) += localK.at(0).at(1);
    K.at(index1).at(index3) += localK.at(0).at(2);
    K.at(index1).at(index4) += localK.at(0).at(3);
    K.at(index1).at(index5) += localK.at(0).at(4);
    K.at(index1).at(index6) += localK.at(0).at(5);
    K.at(index1).at(index7) += localK.at(0).at(6);
    K.at(index1).at(index8) += localK.at(0).at(7);
    K.at(index1).at(index9) += localK.at(0).at(8);

    K.at(index2).at(index1) += localK.at(1).at(0);
    K.at(index2).at(index2) += localK.at(1).at(1);
    K.at(index2).at(index3) += localK.at(1).at(2);
    K.at(index2).at(index4) += localK.at(1).at(3);
    K.at(index2).at(index5) += localK.at(1).at(4);
    K.at(index2).at(index6) += localK.at(1).at(5);
    K.at(index2).at(index7) += localK.at(1).at(6);
    K.at(index2).at(index8) += localK.at(1).at(7);
    K.at(index2).at(index9) += localK.at(1).at(8);

    K.at(index3).at(index1) += localK.at(2).at(0);
    K.at(index3).at(index2) += localK.at(2).at(1);
    K.at(index3).at(index3) += localK.at(2).at(2);
    K.at(index3).at(index4) += localK.at(2).at(3);
    K.at(index3).at(index5) += localK.at(2).at(4);
    K.at(index3).at(index6) += localK.at(2).at(5);
    K.at(index3).at(index7) += localK.at(2).at(6);
    K.at(index3).at(index8) += localK.at(2).at(7);
    K.at(index3).at(index9) += localK.at(2).at(8);

    K.at(index4).at(index1) += localK.at(3).at(0);
    K.at(index4).at(index2) += localK.at(3).at(1);
    K.at(index4).at(index3) += localK.at(3).at(2);
    K.at(index4).at(index4) += localK.at(3).at(3);
    K.at(index4).at(index5) += localK.at(3).at(4);
    K.at(index4).at(index6) += localK.at(3).at(5);
    K.at(index4).at(index7) += localK.at(3).at(6);
    K.at(index4).at(index8) += localK.at(3).at(7);
    K.at(index4).at(index9) += localK.at(3).at(8);

    K.at(index5).at(index1) += localK.at(4).at(0);
    K.at(index5).at(index2) += localK.at(4).at(1);
    K.at(index5).at(index3) += localK.at(4).at(2);
    K.at(index5).at(index4) += localK.at(4).at(3);
    K.at(index5).at(index5) += localK.at(4).at(4);
    K.at(index5).at(index6) += localK.at(4).at(5);
    K.at(index5).at(index7) += localK.at(4).at(6);
    K.at(index5).at(index8) += localK.at(4).at(7);
    K.at(index5).at(index9) += localK.at(4).at(8);

    K.at(index6).at(index1) += localK.at(5).at(0);
    K.at(index6).at(index2) += localK.at(5).at(1);
    K.at(index6).at(index3) += localK.at(5).at(2);
    K.at(index6).at(index4) += localK.at(5).at(3);
    K.at(index6).at(index5) += localK.at(5).at(4);
    K.at(index6).at(index6) += localK.at(5).at(5);
    K.at(index6).at(index7) += localK.at(5).at(6);
    K.at(index6).at(index8) += localK.at(5).at(7);
    K.at(index6).at(index9) += localK.at(5).at(8);

    K.at(index7).at(index1) += localK.at(6).at(0);
    K.at(index7).at(index2) += localK.at(6).at(1);
    K.at(index7).at(index3) += localK.at(6).at(2);
    K.at(index7).at(index4) += localK.at(6).at(3);
    K.at(index7).at(index5) += localK.at(6).at(4);
    K.at(index7).at(index6) += localK.at(6).at(5);
    //K.at(index7).at(index7) += localK.at(6).at(6);
    //K.at(index7).at(index8) += localK.at(6).at(7);
    //K.at(index7).at(index9) += localK.at(6).at(8);

    K.at(index8).at(index1) += localK.at(7).at(0);
    K.at(index8).at(index2) += localK.at(7).at(1);
    K.at(index8).at(index3) += localK.at(7).at(2);
    K.at(index8).at(index4) += localK.at(7).at(3);
    K.at(index8).at(index5) += localK.at(7).at(4);
    K.at(index8).at(index6) += localK.at(7).at(5);
    //K.at(index8).at(index7) += localK.at(7).at(6);
    //K.at(index8).at(index8) += localK.at(7).at(7);
    //K.at(index8).at(index9) += localK.at(7).at(8);

    K.at(index9).at(index1) += localK.at(8).at(0);
    K.at(index9).at(index2) += localK.at(8).at(1);
    K.at(index9).at(index3) += localK.at(8).at(2);
    K.at(index9).at(index4) += localK.at(8).at(3);
    K.at(index9).at(index5) += localK.at(8).at(4);
    K.at(index9).at(index6) += localK.at(8).at(5);
    //K.at(index9).at(index7) += localK.at(8).at(6);
    //K.at(index9).at(index8) += localK.at(8).at(7);
    //K.at(index9).at(index9) += localK.at(8).at(8);
}

void assemblyb(element e,Vector localb,Vector &b,int nnodes){
    int index1 = e.getNode1() - 1;
    int index2 = e.getNode2() - 1;
    int index3 = e.getNode3() - 1;
    int index4 = index1+nnodes, index5 = index2+nnodes, index6 = index3+nnodes;

    b.at(index1) += localb.at(0);
    b.at(index2) += localb.at(1);
    b.at(index3) += localb.at(2);
    b.at(index4) += localb.at(3);
    b.at(index5) += localb.at(4);
    b.at(index6) += localb.at(5);
}

