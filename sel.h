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
    int index4 = e.getNode4() - 1;

    int index5 = index1 + nnodes;
    int index6 = index2 + nnodes;
    int index7 = index3 + nnodes;
    int index8 = index4 + nnodes;

    int index9  = index1 + 2*nnodes;
    int index10 = index2 + 2*nnodes;
    int index11 = index3 + 2*nnodes;
    int index12 = index4 + 2*nnodes;

    int index13 = index1 + 3*nnodes;
    int index14 = index2 + 3*nnodes;
    int index15 = index3 + 3*nnodes;
    int index16 = index4 + 3*nnodes;

    K.at(index1).at(index1)  += localK.at(0).at(0);
    K.at(index1).at(index2)  += localK.at(0).at(1);
    K.at(index1).at(index3)  += localK.at(0).at(2);
    K.at(index1).at(index4)  += localK.at(0).at(3);
    K.at(index1).at(index5)  += localK.at(0).at(4);
    K.at(index1).at(index6)  += localK.at(0).at(5);
    K.at(index1).at(index7)  += localK.at(0).at(6);
    K.at(index1).at(index8)  += localK.at(0).at(7);
    K.at(index1).at(index9)  += localK.at(0).at(8);
    K.at(index1).at(index10) += localK.at(0).at(9);
    K.at(index1).at(index11) += localK.at(0).at(10);
    K.at(index1).at(index12) += localK.at(0).at(11);
    K.at(index1).at(index13) += localK.at(0).at(12);
    K.at(index1).at(index14) += localK.at(0).at(13);
    K.at(index1).at(index15) += localK.at(0).at(14);
    K.at(index1).at(index16) += localK.at(0).at(15);

    K.at(index2).at(index1)  += localK.at(1).at(0);
    K.at(index2).at(index2)  += localK.at(1).at(1);
    K.at(index2).at(index3)  += localK.at(1).at(2);
    K.at(index2).at(index4)  += localK.at(1).at(3);
    K.at(index2).at(index5)  += localK.at(1).at(4);
    K.at(index2).at(index6)  += localK.at(1).at(5);
    K.at(index2).at(index7)  += localK.at(1).at(6);
    K.at(index2).at(index8)  += localK.at(1).at(7);
    K.at(index2).at(index9)  += localK.at(1).at(8);
    K.at(index2).at(index10) += localK.at(1).at(9);
    K.at(index2).at(index11) += localK.at(1).at(10);
    K.at(index2).at(index12) += localK.at(1).at(11);
    K.at(index2).at(index13) += localK.at(1).at(12);
    K.at(index2).at(index14) += localK.at(1).at(13);
    K.at(index2).at(index15) += localK.at(1).at(14);
    K.at(index2).at(index16) += localK.at(1).at(15);

    K.at(index3).at(index1)  += localK.at(2).at(0);
    K.at(index3).at(index2)  += localK.at(2).at(1);
    K.at(index3).at(index3)  += localK.at(2).at(2);
    K.at(index3).at(index4)  += localK.at(2).at(3);
    K.at(index3).at(index5)  += localK.at(2).at(4);
    K.at(index3).at(index6)  += localK.at(2).at(5);
    K.at(index3).at(index7)  += localK.at(2).at(6);
    K.at(index3).at(index8)  += localK.at(2).at(7);
    K.at(index3).at(index9)  += localK.at(2).at(8);
    K.at(index3).at(index10) += localK.at(2).at(9);
    K.at(index3).at(index11) += localK.at(2).at(10);
    K.at(index3).at(index12) += localK.at(2).at(11);
    K.at(index3).at(index13) += localK.at(2).at(12);
    K.at(index3).at(index14) += localK.at(2).at(13);
    K.at(index3).at(index15) += localK.at(2).at(14);
    K.at(index3).at(index16) += localK.at(2).at(15);

    K.at(index4).at(index1)  += localK.at(3).at(0);
    K.at(index4).at(index2)  += localK.at(3).at(1);
    K.at(index4).at(index3)  += localK.at(3).at(2);
    K.at(index4).at(index4)  += localK.at(3).at(3);
    K.at(index4).at(index5)  += localK.at(3).at(4);
    K.at(index4).at(index6)  += localK.at(3).at(5);
    K.at(index4).at(index7)  += localK.at(3).at(6);
    K.at(index4).at(index8)  += localK.at(3).at(7);
    K.at(index4).at(index9)  += localK.at(3).at(8);
    K.at(index4).at(index10) += localK.at(3).at(9);
    K.at(index4).at(index11) += localK.at(3).at(10);
    K.at(index4).at(index12) += localK.at(3).at(11);
    K.at(index4).at(index13) += localK.at(3).at(12);
    K.at(index4).at(index14) += localK.at(3).at(13);
    K.at(index4).at(index15) += localK.at(3).at(14);
    K.at(index4).at(index16) += localK.at(3).at(15);

    K.at(index5).at(index1)  += localK.at(4).at(0);
    K.at(index5).at(index2)  += localK.at(4).at(1);
    K.at(index5).at(index3)  += localK.at(4).at(2);
    K.at(index5).at(index4)  += localK.at(4).at(3);
    K.at(index5).at(index5)  += localK.at(4).at(4);
    K.at(index5).at(index6)  += localK.at(4).at(5);
    K.at(index5).at(index7)  += localK.at(4).at(6);
    K.at(index5).at(index8)  += localK.at(4).at(7);
    K.at(index5).at(index9)  += localK.at(4).at(8);
    K.at(index5).at(index10) += localK.at(4).at(9);
    K.at(index5).at(index11) += localK.at(4).at(10);
    K.at(index5).at(index12) += localK.at(4).at(11);
    K.at(index5).at(index13) += localK.at(4).at(12);
    K.at(index5).at(index14) += localK.at(4).at(13);
    K.at(index5).at(index15) += localK.at(4).at(14);
    K.at(index5).at(index16) += localK.at(4).at(15);

    K.at(index6).at(index1)  += localK.at(5).at(0);
    K.at(index6).at(index2)  += localK.at(5).at(1);
    K.at(index6).at(index3)  += localK.at(5).at(2);
    K.at(index6).at(index4)  += localK.at(5).at(3);
    K.at(index6).at(index5)  += localK.at(5).at(4);
    K.at(index6).at(index6)  += localK.at(5).at(5);
    K.at(index6).at(index7)  += localK.at(5).at(6);
    K.at(index6).at(index8)  += localK.at(5).at(7);
    K.at(index6).at(index9)  += localK.at(5).at(8);
    K.at(index6).at(index10) += localK.at(5).at(9);
    K.at(index6).at(index11) += localK.at(5).at(10);
    K.at(index6).at(index12) += localK.at(5).at(11);
    K.at(index6).at(index13) += localK.at(5).at(12);
    K.at(index6).at(index14) += localK.at(5).at(13);
    K.at(index6).at(index15) += localK.at(5).at(14);
    K.at(index6).at(index16) += localK.at(5).at(15);

    K.at(index7).at(index1)  += localK.at(6).at(0);
    K.at(index7).at(index2)  += localK.at(6).at(1);
    K.at(index7).at(index3)  += localK.at(6).at(2);
    K.at(index7).at(index4)  += localK.at(6).at(3);
    K.at(index7).at(index5)  += localK.at(6).at(4);
    K.at(index7).at(index6)  += localK.at(6).at(5);
    K.at(index7).at(index7)  += localK.at(6).at(6);
    K.at(index7).at(index8)  += localK.at(6).at(7);
    K.at(index7).at(index9)  += localK.at(6).at(8);
    K.at(index7).at(index10) += localK.at(6).at(9);
    K.at(index7).at(index11) += localK.at(6).at(10);
    K.at(index7).at(index12) += localK.at(6).at(11);
    K.at(index7).at(index13) += localK.at(6).at(12);
    K.at(index7).at(index14) += localK.at(6).at(13);
    K.at(index7).at(index15) += localK.at(6).at(14);
    K.at(index7).at(index16) += localK.at(6).at(15);

    K.at(index8).at(index1)  += localK.at(7).at(0);
    K.at(index8).at(index2)  += localK.at(7).at(1);
    K.at(index8).at(index3)  += localK.at(7).at(2);
    K.at(index8).at(index4)  += localK.at(7).at(3);
    K.at(index8).at(index5)  += localK.at(7).at(4);
    K.at(index8).at(index6)  += localK.at(7).at(5);
    K.at(index8).at(index7)  += localK.at(7).at(6);
    K.at(index8).at(index8)  += localK.at(7).at(7);
    K.at(index8).at(index9)  += localK.at(7).at(8);
    K.at(index8).at(index10) += localK.at(7).at(9);
    K.at(index8).at(index11) += localK.at(7).at(10);
    K.at(index8).at(index12) += localK.at(7).at(11);
    K.at(index8).at(index13) += localK.at(7).at(12);
    K.at(index8).at(index14) += localK.at(7).at(13);
    K.at(index8).at(index15) += localK.at(7).at(14);
    K.at(index8).at(index16) += localK.at(7).at(15);

    K.at(index9).at(index1)  += localK.at(8).at(0);
    K.at(index9).at(index2)  += localK.at(8).at(1);
    K.at(index9).at(index3)  += localK.at(8).at(2);
    K.at(index9).at(index4)  += localK.at(8).at(3);
    K.at(index9).at(index5)  += localK.at(8).at(4);
    K.at(index9).at(index6)  += localK.at(8).at(5);
    K.at(index9).at(index7)  += localK.at(8).at(6);
    K.at(index9).at(index8)  += localK.at(8).at(7);
    K.at(index9).at(index9)  += localK.at(8).at(8);
    K.at(index9).at(index10) += localK.at(8).at(9);
    K.at(index9).at(index11) += localK.at(8).at(10);
    K.at(index9).at(index12) += localK.at(8).at(11);
    K.at(index9).at(index13) += localK.at(8).at(12);
    K.at(index9).at(index14) += localK.at(8).at(13);
    K.at(index9).at(index15) += localK.at(8).at(14);
    K.at(index9).at(index16) += localK.at(8).at(15);

    K.at(index10).at(index1)  += localK.at(9).at(0);
    K.at(index10).at(index2)  += localK.at(9).at(1);
    K.at(index10).at(index3)  += localK.at(9).at(2);
    K.at(index10).at(index4)  += localK.at(9).at(3);
    K.at(index10).at(index5)  += localK.at(9).at(4);
    K.at(index10).at(index6)  += localK.at(9).at(5);
    K.at(index10).at(index7)  += localK.at(9).at(6);
    K.at(index10).at(index8)  += localK.at(9).at(7);
    K.at(index10).at(index9)  += localK.at(9).at(8);
    K.at(index10).at(index10) += localK.at(9).at(9);
    K.at(index10).at(index11) += localK.at(9).at(10);
    K.at(index10).at(index12) += localK.at(9).at(11);
    K.at(index10).at(index13) += localK.at(9).at(12);
    K.at(index10).at(index14) += localK.at(9).at(13);
    K.at(index10).at(index15) += localK.at(9).at(14);
    K.at(index10).at(index16) += localK.at(9).at(15);

    K.at(index11).at(index1)  += localK.at(10).at(0);
    K.at(index11).at(index2)  += localK.at(10).at(1);
    K.at(index11).at(index3)  += localK.at(10).at(2);
    K.at(index11).at(index4)  += localK.at(10).at(3);
    K.at(index11).at(index5)  += localK.at(10).at(4);
    K.at(index11).at(index6)  += localK.at(10).at(5);
    K.at(index11).at(index7)  += localK.at(10).at(6);
    K.at(index11).at(index8)  += localK.at(10).at(7);
    K.at(index11).at(index9)  += localK.at(10).at(8);
    K.at(index11).at(index10) += localK.at(10).at(9);
    K.at(index11).at(index11) += localK.at(10).at(10);
    K.at(index11).at(index12) += localK.at(10).at(11);
    K.at(index11).at(index13) += localK.at(10).at(12);
    K.at(index11).at(index14) += localK.at(10).at(13);
    K.at(index11).at(index15) += localK.at(10).at(14);
    K.at(index11).at(index16) += localK.at(10).at(15);

    K.at(index12).at(index1)  += localK.at(11).at(0);
    K.at(index12).at(index2)  += localK.at(11).at(1);
    K.at(index12).at(index3)  += localK.at(11).at(2);
    K.at(index12).at(index4)  += localK.at(11).at(3);
    K.at(index12).at(index5)  += localK.at(11).at(4);
    K.at(index12).at(index6)  += localK.at(11).at(5);
    K.at(index12).at(index7)  += localK.at(11).at(6);
    K.at(index12).at(index8)  += localK.at(11).at(7);
    K.at(index12).at(index9)  += localK.at(11).at(8);
    K.at(index12).at(index10) += localK.at(11).at(9);
    K.at(index12).at(index11) += localK.at(11).at(10);
    K.at(index12).at(index12) += localK.at(11).at(11);
    K.at(index12).at(index13) += localK.at(11).at(12);
    K.at(index12).at(index14) += localK.at(11).at(13);
    K.at(index12).at(index15) += localK.at(11).at(14);
    K.at(index12).at(index16) += localK.at(11).at(15);

    K.at(index13).at(index1)  += localK.at(12).at(0);
    K.at(index13).at(index2)  += localK.at(12).at(1);
    K.at(index13).at(index3)  += localK.at(12).at(2);
    K.at(index13).at(index4)  += localK.at(12).at(3);
    K.at(index13).at(index5)  += localK.at(12).at(4);
    K.at(index13).at(index6)  += localK.at(12).at(5);
    K.at(index13).at(index7)  += localK.at(12).at(6);
    K.at(index13).at(index8)  += localK.at(12).at(7);
    K.at(index13).at(index9)  += localK.at(12).at(8);
    K.at(index13).at(index10) += localK.at(12).at(9);
    K.at(index13).at(index11) += localK.at(12).at(10);
    K.at(index13).at(index12) += localK.at(12).at(11);
    K.at(index13).at(index13) += localK.at(12).at(12);
    K.at(index13).at(index14) += localK.at(12).at(13);
    K.at(index13).at(index15) += localK.at(12).at(14);
    K.at(index13).at(index16) += localK.at(12).at(15);

    K.at(index14).at(index1)  += localK.at(13).at(0);
    K.at(index14).at(index2)  += localK.at(13).at(1);
    K.at(index14).at(index3)  += localK.at(13).at(2);
    K.at(index14).at(index4)  += localK.at(13).at(3);
    K.at(index14).at(index5)  += localK.at(13).at(4);
    K.at(index14).at(index6)  += localK.at(13).at(5);
    K.at(index14).at(index7)  += localK.at(13).at(6);
    K.at(index14).at(index8)  += localK.at(13).at(7);
    K.at(index14).at(index9)  += localK.at(13).at(8);
    K.at(index14).at(index10) += localK.at(13).at(9);
    K.at(index14).at(index11) += localK.at(13).at(10);
    K.at(index14).at(index12) += localK.at(13).at(11);
    K.at(index14).at(index13) += localK.at(13).at(12);
    K.at(index14).at(index14) += localK.at(13).at(13);
    K.at(index14).at(index15) += localK.at(13).at(14);
    K.at(index14).at(index16) += localK.at(13).at(15);

    K.at(index15).at(index1)  += localK.at(14).at(0);
    K.at(index15).at(index2)  += localK.at(14).at(1);
    K.at(index15).at(index3)  += localK.at(14).at(2);
    K.at(index15).at(index4)  += localK.at(14).at(3);
    K.at(index15).at(index5)  += localK.at(14).at(4);
    K.at(index15).at(index6)  += localK.at(14).at(5);
    K.at(index15).at(index7)  += localK.at(14).at(6);
    K.at(index15).at(index8)  += localK.at(14).at(7);
    K.at(index15).at(index9)  += localK.at(14).at(8);
    K.at(index15).at(index10) += localK.at(14).at(9);
    K.at(index15).at(index11) += localK.at(14).at(10);
    K.at(index15).at(index12) += localK.at(14).at(11);
    K.at(index15).at(index13) += localK.at(14).at(12);
    K.at(index15).at(index14) += localK.at(14).at(13);
    K.at(index15).at(index15) += localK.at(14).at(14);
    K.at(index15).at(index16) += localK.at(14).at(15);

    K.at(index16).at(index1)  += localK.at(15).at(0);
    K.at(index16).at(index2)  += localK.at(15).at(1);
    K.at(index16).at(index3)  += localK.at(15).at(2);
    K.at(index16).at(index4)  += localK.at(15).at(3);
    K.at(index16).at(index5)  += localK.at(15).at(4);
    K.at(index16).at(index6)  += localK.at(15).at(5);
    K.at(index16).at(index7)  += localK.at(15).at(6);
    K.at(index16).at(index8)  += localK.at(15).at(7);
    K.at(index16).at(index9)  += localK.at(15).at(8);
    K.at(index16).at(index10) += localK.at(15).at(9);
    K.at(index16).at(index11) += localK.at(15).at(10);
    K.at(index16).at(index12) += localK.at(15).at(11);
    K.at(index16).at(index13) += localK.at(15).at(12);
    K.at(index16).at(index14) += localK.at(15).at(13);
    K.at(index16).at(index15) += localK.at(15).at(14);
    K.at(index16).at(index16) += localK.at(15).at(15);

}

void assemblyb(element e,Vector localb,Vector &b,int nnodes){
    int index1 = e.getNode1() - 1;
    int index2 = e.getNode2() - 1;
    int index3 = e.getNode3() - 1;
    int index4 = e.getNode4() - 1;

    int index5 = index1 + nnodes;
    int index6 = index2 + nnodes;
    int index7 = index3 + nnodes;
    int index8 = index4 + nnodes;

    int index9  = index1 + 2*nnodes;
    int index10 = index2 + 2*nnodes;
    int index11 = index3 + 2*nnodes;
    int index12 = index4 + 2*nnodes;    

    int index13 = index1 + 3*nnodes;
    int index14 = index2 + 3*nnodes;
    int index15 = index3 + 3*nnodes;
    int index16 = index4 + 3*nnodes;    


    b.at(index1) += localb.at(0);
    b.at(index2) += localb.at(1);
    b.at(index3) += localb.at(2);
    b.at(index4) += localb.at(3);
    b.at(index5) += localb.at(4);
    b.at(index6) += localb.at(5);
    b.at(index7) += localb.at(6);
    b.at(index8) += localb.at(7);
    b.at(index9) += localb.at(8);
    b.at(index10) += localb.at(9);
    b.at(index11) += localb.at(10);
    b.at(index12) += localb.at(11);
    b.at(index13) += localb.at(12);
    b.at(index14) += localb.at(13);
    b.at(index15) += localb.at(14);
    b.at(index16) += localb.at(15);

}

