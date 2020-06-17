void crearSistemasLocales(mesh &m,vector<Matrix> &localKs,vector<Vector> &localbs){
    for(int i=0;i<m.getSize(ELEMENTS);i++){
        localKs.push_back(createLocalM(i,m));
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

