//======================================================================================================
// File pySOQCS.cpp
//
// INTERFACE WITH PYTHON LIBRARY. C++ SIDE.
//
// Copyright © 2022 National University of Ireland Maynooth, Maynooth University. All rights reserved.
// The contents of this file are subject to the licence terms detailed in LICENCE.TXT available in the
// root directory of this source tree. Use of the source code in this file is only permitted under the
// terms of the licence in LICENCE.TXT.
//======================================================================================================

#include "./src/soqcs.h"


//----------------------------------------
//
//  Free memory of arrays
//
//----------------------------------------
void free_mem(char *mem){
    delete mem;
}


//----------------------------------------
//
//  Converter of array to matrix of integers
//
//----------------------------------------

mati to_mati(int *int_array, int n, int m){
    int i;
    int j;
    mati cnv;

    cnv.resize(n,m);
    for(i=0;i<n;i++){
        for(j=0;j<m;j++){
           cnv(i,j)=int_array[i*m+j];
        }
    }
    return cnv;
}


//----------------------------------------
//
//  Converter of array to matrix of doubles
//
//----------------------------------------
matd to_matrix(double *double_array, int n, int m){
    int i;
    int j;
    matd cnv;

    cnv.resize(n,m);
    for(i=0;i<n;i++){
        for(j=0;j<m;j++){
           cnv(i,j)=double_array[i*m+j];
        }
    }
    return cnv;
}


//----------------------------------------
//
//  Converter of array to vector of integers
//
//----------------------------------------
veci to_veci(int *int_array, int n){
    int i;
    veci cnv;

    cnv.resize(n);
    for(i=0;i<n;i++){
        cnv(i)=int_array[i];
    }
    return cnv;
}


//----------------------------------------
//
//  Converter of array to matrix of doubles
//
//----------------------------------------
double *to_dptr(matc mtx){
    int i;
    int j;
    int k;
    double *cnv;

    cnv=new double[2*mtx.rows()*mtx.cols()];
//    cout << "Rows:" << mtx.rows() << endl;
//    cout << "Cols:" << mtx.cols() << endl;

    k=0;
    for(i=0;i<mtx.rows();i++){
        for(j=0;j<mtx.cols();j++){
            cnv[k] = real(mtx(i,j));
            cnv[k+1]=imag(mtx(i,j));
            k=k+2;
//            cout << "k:" << k << endl;
        }
    }
    return cnv;
}


//----------------------------------------
//
//  Interface with Python. C++/C side.
//
//----------------------------------------
extern "C" {
    //--------------------------------------------------------------------------------------------------------------------------
    // GENERAL
    // Interface support
    void free_ptr(char *mem){ free_mem(mem); }
    // Configure SOQCS
    void all_cfg_soqcs(int nph){cfg_soqcs(nph);}
    //--------------------------------------------------------------------------------------------------------------------------


    //--------------------------------------------------------------------------------------------------------------------------
    // QOCIRCUIT
    // Management functions
    long int qoc_new_qocircuit(int i_nch, int i_nm, int i_ns, int i_np,  double i_dtp, int clock, int i_R, bool loss, int ikind){ char ckind='G';if(ikind==1) ckind='E';
                                                                                                  return (long int)new qocircuit(i_nch,i_nm,i_ns, i_np, i_dtp, clock,i_R,loss,ckind);}
    void qoc_destroy_qocircuit(long int qoc){qocircuit* aux=(qocircuit*)qoc; delete aux; }
    int qoc_num_levels(long int qoc){qocircuit* aux=(qocircuit*)qoc; return aux->num_levels();}

    // Circuit elements
    //      Basic elements
    void qoc_random_circuit(long int qoc){ qocircuit* aux=(qocircuit*)qoc; aux->random_circuit();}
    void qoc_NSX(long int  qoc, int i_ch1, int i_ch2, int i_ch3){ qocircuit* aux=(qocircuit*)qoc; aux->NSX(i_ch1,i_ch2,i_ch3);}
    void qoc_beamsplitter(long int qoc, int i_ch1, int i_ch2, double theta, double phi){qocircuit* aux=(qocircuit*)qoc; aux->beamsplitter(i_ch1,i_ch2, theta, phi);}
    void qoc_dielectric(long int qoc, int i_ch1, int i_ch2, double ret, double imt, double rer, double imr){ cmplx t=ret+jm*imt; cmplx r=rer+jm*imr; qocircuit* aux=(qocircuit*)qoc; aux->dielectric(i_ch1,i_ch2,t,r);}
    void qoc_MMI2(long int qoc, int i_ch1, int i_ch2){ qocircuit* aux=(qocircuit*)qoc; aux->MMI2(i_ch1, i_ch2);}
    void qoc_rewire(long int qoc,int i_ch1,int i_ch2){ qocircuit* aux=(qocircuit*)qoc; aux->rewire(i_ch1, i_ch2);}
    void qoc_phase_shifter(long int qoc, int i_ch, double ret, double imt){ cmplx t= ret+jm*imt; qocircuit* aux=(qocircuit*)qoc; aux->phase_shifter(i_ch,t);}

    // Polarization elements
    void qoc_rotator(long int qoc, int i_ch, double theta, double phi){ qocircuit* aux=(qocircuit*)qoc; aux->rotator(i_ch,theta,phi);}
    void qoc_pol_beamsplitter(long int qoc, int i_ch1, int i_ch2, int P){qocircuit* aux=(qocircuit*)qoc; aux->pol_beamsplitter(i_ch1,i_ch2, P);}
    void qoc_pol_phase_shifter(long int qoc, int i_ch, int P, double phi){qocircuit* aux=(qocircuit*)qoc; aux->pol_phase_shifter(i_ch,P, phi);}
    void qoc_pol_filter(long int qoc, int i_ch, int P){qocircuit* aux=(qocircuit*)qoc; aux->pol_filter(i_ch,P);}
    void qoc_waveplate(long int qoc, int i_ch, double alpha, double gamma){qocircuit* aux=(qocircuit*)qoc; aux->waveplate(i_ch, alpha, gamma);}

    //      Detection elements
    void qoc_detector(long int  qoc, int i_ch, int cond, int pol, int mpi, int mpo, double eff, double blnk, double gamma){ qocircuit* aux=(qocircuit*)qoc; aux->detector(i_ch,cond,pol,mpi,mpo,eff,blnk,gamma);}
    void qoc_noise(long int qoc, double stdev2){ qocircuit* aux=(qocircuit*)qoc; aux->noise(stdev2);}

    //      Emitter and distinguishability model.
    int qoc_def_packet(long int qoc, int n, double t, double f, double w){ qocircuit* aux=(qocircuit*)qoc;
                                                                           return aux->def_packet(n,t, f, w);}
    double qoc_emitted_vis(long int qoc, int i,int j){ qocircuit* aux=(qocircuit*)qoc; return aux->emitted_vis(i, j);}
    void   qoc_emitter(long int qoc){ qocircuit* aux=(qocircuit*)qoc; aux->emitter ();}
    void   qoc_delay(long int qoc, int ch){ qocircuit* aux=(qocircuit*)qoc; aux->delay(ch);}

    // Print functions
    void qoc_prnt(long int qoc, int format){ qocircuit* aux=(qocircuit*)qoc; aux->prnt(format); cout << flush;}
    //--------------------------------------------------------------------------------------------------------------------------


    //--------------------------------------------------------------------------------------------------------------------------
    // STATE
    // Management functions
    long int st_new_state(int i_level, int i_maxket){ return (long int)new state(i_level,i_maxket);}
    void st_destroy_state(long int st){state* aux=(state*)st; delete aux; }

    // State manipulation methods.
    double *st_braket(long int st1,long int st2){ state *auxst1=(state*)st1;
                                                  state *auxst2=(state*)st2;
                                                  cmplx value=auxst1->braket(auxst2);
                                                  matc matvalue;
                                                  matvalue.resize(1,2);
                                                  matvalue(0,0)=real(value);
                                                  matvalue(0,1)=imag(value);
                                                  double *auxdouble=to_dptr(matvalue);
                                                  return auxdouble;}

    void st_normalize(long int st){ state *auxst=(state*)st; auxst->normalize(); }
    void st_rephase(long int st, int *term, int n, int m, long int qoc){ state *auxst=(state*)st;
                                                                         qocircuit *auxqoc=(qocircuit*)qoc;
                                                                         mati imat=to_mati(term,n,m);
                                                                         auxst->rephase(imat,auxqoc); }

    void st_add_raw_term(long int st, double rampl, double iampl, int *term){ state *auxst=(state*)st;
                                                                              cmplx ampl=rampl+jm*iampl;
                                                                              auxst->add_term(ampl,term); }

    void st_add_term(long int st, double rampl, double iampl, int *term, int n, int m, long int qoc){ state *auxst=(state*)st;
                                                                                                      qocircuit *auxqoc=(qocircuit*)qoc;
                                                                                                      cmplx ampl=rampl+jm*iampl;
                                                                                                      mati imat=to_mati(term,n,m);
                                                                                                      auxst->add_term(ampl,imat,auxqoc); }

    long int st_post_selection(long int st, long int prj){ state* auxst=(state *)st; projector* auxprj=(projector*)prj; return (long int)auxst->post_selection(auxprj); }

    // Print methods
    void  st_prnt_state(long int st, int format, int column, bool loss, long int qoc){  state* auxst=(state*)st; qocircuit* auxqoc=(qocircuit*)qoc; auxst->prnt_state(format,column,loss,auxqoc); cout << flush;}

    // Qubit codification methods
    long int st_encode(long int st, int *qdef, int nqbits, long int qoc)     {state *auxst=(state *)st;
                                                                              qocircuit *auxqoc=(qocircuit*)qoc;
                                                                              mati imat=to_mati(qdef,2,nqbits);
                                                                              state*auxencoded=auxst->encode(imat,auxqoc);
                                                                              return (long int) auxencoded;}

    long int st_decode(long int st, int *qdef, int nqbits, long int ancilla, long int qoc)     {state *auxst=(state *)st;
                                                                                                state *auxan=(state *)ancilla;
                                                                                                qocircuit *auxqoc=(qocircuit*)qoc;
                                                                                                mati imat=to_mati(qdef,2,nqbits);
                                                                                                state*auxdecoded=auxst->decode(imat,auxan,auxqoc);
                                                                                                return (long int) auxdecoded;}

    long int st_pol_encode(long int st, int *qdef, int nqbits, long int qoc)     {state *auxst=(state *)st;
                                                                                  qocircuit *auxqoc=(qocircuit*)qoc;
                                                                                  veci ivec=to_veci(qdef,nqbits);
                                                                                  state*auxencoded=auxst->pol_encode(ivec,auxqoc);
                                                                                  return (long int) auxencoded;}

    long int st_pol_decode(long int st, int *qdef, int nqbits, long int ancilla, long int qoc)     {state *auxst=(state *)st;
                                                                                                    state *auxan=(state *)ancilla;
                                                                                                    qocircuit *auxqoc=(qocircuit*)qoc;
                                                                                                    veci ivec=to_veci(qdef,nqbits);
                                                                                                    state*auxdecoded=auxst->pol_decode(ivec,auxan,auxqoc);
                                                                                                    return (long int) auxdecoded;}

    // PROJECTOR
    long int prj_new_projector(int i_level, int i_maxket){ return (long int)new projector(i_level,i_maxket);}
    void prj_destroy_projector(long int prj){projector* aux=(projector *)prj; delete aux; }

    // State manipulation methods.
    void prj_add_term(long int prj, double rampl, double iampl, int *term, int n, int m, long int qoc){ projector *auxprj=(projector *)prj;
                                                                                                        qocircuit *auxqoc=(qocircuit*)qoc;
                                                                                                        cmplx ampl=rampl+jm*iampl;
                                                                                                        mati imat;
                                                                                                        imat=to_mati(term,n,m);
                                                                                                        auxprj->add_term(ampl,imat,auxqoc);}
    //--------------------------------------------------------------------------------------------------------------------------


    //--------------------------------------------------------------------------------------------------------------------------
    // PBINS
    // Management functions
    long int pb_new_pbin(int i_level, int i_maxket){ return (long int)new p_bin(i_level,i_maxket);}
    void pb_destroy_pbin(long int pbin){p_bin* aux=(p_bin *)pbin; delete aux; }

    // Update pdate methods
    void pb_add_state(long int pbin, long int st){p_bin *auxpb=(p_bin*)pbin; state *auxst=(state*)st; auxpb->add_state(auxst);}
    double pb_trace(long int pbin){p_bin *auxpb=(p_bin*)pbin; return auxpb->trace();}
    void pb_normalize(long int pbin){p_bin *auxpb=(p_bin*)pbin; auxpb->normalize();}

    // Manipulation methods
    long int pb_calc_measure(long int pbin, long int qoc){p_bin *auxpb=(p_bin *)pbin; qocircuit *auxqoc=(qocircuit *)qoc; return (long int)auxpb->calc_measure(auxqoc);}

    //Bin consultation methods
    int pb_nbins(long int pbin){p_bin *auxpb=(p_bin*)pbin; return auxpb->nket;}
    int pb_num_levels(long int pbin){p_bin *auxpb=(p_bin*)pbin; return auxpb->nlevel;}
    char *pb_tag(long int pbin, int index){p_bin *auxpb=(p_bin*)pbin; string value =auxpb->tag(index); char* char_array=new char[value.length()+1]; strcpy(char_array, value.c_str()); return char_array;}
    double pb_prob(long int pbin,  int index){p_bin *auxpb=(p_bin*)pbin; return auxpb->prob(index);}
    double pb_prob_def_qoc(long int pbin, int *def,int n, int m, long int qoc){  p_bin *auxpb=(p_bin *) pbin;
                                                                                 qocircuit *auxqoc=(qocircuit*)qoc;
                                                                                 mati imat;
                                                                                 imat=to_mati(def,n,m);
                                                                                 return auxpb->prob(imat,auxqoc);
                                                                              }

    double pb_prob_def(long int pbin, int *def,int n, int m, long int dev){      p_bin *auxpb=(p_bin *) pbin;
                                                                                 qodev *auxdev=(qodev*)dev;
                                                                                 mati imat;
                                                                                 imat=to_mati(def,n,m);
                                                                                 return auxpb->prob(imat,auxdev);
                                                                              }

    // Print bins
    void  pb_prnt_bins_qoc(long int pbin, int format, double thresh, bool loss, long int qoc){p_bin *auxpb=(p_bin*)pbin; qocircuit *auxqoc=(qocircuit*)qoc; auxpb->prnt_bins(format,thresh,loss,auxqoc); cout << flush;}
    void  pb_prnt_bins(long int pbin, int format, double thresh, bool loss, long int dev){p_bin *auxpb=(p_bin*)pbin; qodev *auxdev=(qodev*)dev; auxpb->prnt_bins(format,thresh,loss,auxdev); cout << flush;}

    // Qubit codification methods
    long int pb_translate(long int pbin, int *qdef, int nqbits, long int dev){p_bin *auxpb=(p_bin *)pbin;
                                                                              qodev *auxdev=(qodev *)dev;
                                                                              mati imat=to_mati(qdef,2,nqbits);
                                                                              return (long int)auxpb->translate(imat, auxdev->circ);}

    long int pb_pol_translate(long int pbin, int *qdef, int nqbits, long int dev){p_bin *auxpb=(p_bin *)pbin;
                                                                              qodev *auxdev=(qodev *)dev;
                                                                              veci ivec=to_veci(qdef,nqbits);
                                                                              return (long int)auxpb->pol_translate(ivec, auxdev->circ);}
    //--------------------------------------------------------------------------------------------------------------------------


    //--------------------------------------------------------------------------------------------------------------------------
    // DMAT
    // Management functions
    long int dm_new_dmat(int i_mem){ return (long int)new dmatrix(i_mem);}
    void dm_destroy_dmat(long int dmat){dmatrix* aux=(dmatrix *)dmat; delete aux; }

    // Basic matrix operations
    double dm_trace(long int dmat){dmatrix *auxdm=(dmatrix *) dmat; return auxdm->trace();}
    void dm_normalize(long int dmat){dmatrix *auxdm=(dmatrix *) dmat; auxdm->normalize();}
    double dm_fidelity(long int dmat, long int st){dmatrix *auxdm=(dmatrix *) dmat; state *auxst=(state*)st; return auxdm->fidelity(auxst);}


    // Update pdate methods
    void dm_add_state_qoc(long int dmat, long int st, long int qoc){dmatrix *auxdm=(dmatrix*)dmat; state *auxst=(state*)st;  qocircuit *auxqoc=(qocircuit *)qoc; auxdm->add_state(auxst,auxqoc);}
    void dm_add_state(long int dmat, long int st, long int dev){dmatrix *auxdm=(dmatrix*)dmat; state *auxst=(state*)st;  qodev *auxdev=(qodev *)dev; auxdm->add_state(auxst,auxdev);}
    long int dm_calc_measure_qoc(long int dmat, long int qoc){dmatrix *auxdm=(dmatrix *)dmat; qocircuit *auxqoc=(qocircuit *)qoc; return (long int)auxdm->calc_measure(auxqoc);}
    long int dm_calc_measure(long int dmat, long int dev){dmatrix *auxdm=(dmatrix *)dmat; qodev *auxdev=(qodev *)dev; return (long int)auxdm->calc_measure(auxdev);}

    // Print bins
    void  dm_prnt_mtx_qoc(long int dmat, int format, double thresh, long int qoc){dmatrix *auxdm=(dmatrix*)dmat; qocircuit *auxqoc=(qocircuit*)qoc; auxdm->prnt_mtx(format,thresh,auxqoc); cout << flush;}
    void  dm_prnt_mtx(long int dmat, int format, double thresh, long int dev){dmatrix *auxdm=(dmatrix*)dmat; qodev *auxdev=(qodev *)dev; auxdm->prnt_mtx(format,thresh,auxdev); cout << flush;}


    // Qubit codification methods
    long int dm_translate(long int dmat, int *qdef, int nqbits, long int dev){dmatrix *auxdm=(dmatrix *)dmat;
                                                                              qodev *auxdev=(qodev *)dev;
                                                                              mati imat=to_mati(qdef,2,nqbits);
                                                                              return (long int)auxdm->translate(imat, auxdev->circ);}

    long int dm_pol_translate(long int dmat, int *qdef, int nqbits, long int dev){dmatrix *auxdm=(dmatrix *)dmat;
                                                                              qodev *auxdev=(qodev *)dev;
                                                                              veci ivec=to_veci(qdef,nqbits);
                                                                              return (long int)auxdm->pol_translate(ivec, auxdev->circ);}

    //--------------------------------------------------------------------------------------------------------------------------
    // QODEV
    // Management functions
    long int dev_new_qodev(int i_nph, int i_nch, int i_nm, int i_ns, int i_np, double i_dtp, int clock, int i_R, bool loss, int ikind, int i_maxket){ char ckind='G'; if(ikind==1) ckind='E';
                                                                                  return (long int)new qodev(i_nph,i_nch,i_nm,i_ns,i_np,i_dtp,clock,i_R, loss, ckind, i_maxket);}
    void dev_destroy_qodev(long int dev){qodev *aux=(qodev *)dev; delete aux; }
    void dev_concatenate(long int dev1, long int dev2){qodev *auxdev1=(qodev *)dev1; qodev *auxdev2=(qodev *)dev2; auxdev1->concatenate(auxdev2);}
    void dev_add_gate(long int dev1, int *chlist, int n, long int dev2){qodev *auxdev1=(qodev *)dev1;
                                                                        qodev *auxdev2=(qodev *)dev2;
                                                                        veci ivec=to_veci(chlist,n);
                                                                        auxdev1->add_gate(ivec,auxdev2);
                                                                        }


    // Initial state definition
    int dev_add_photons(long int dev, int N, int ch, int P, double t, double f, double w){  qodev *auxdev=(qodev *) dev; return auxdev->add_photons(N,ch,P,t,f,w);}
    void dev_add_QD(long int dev, int ch1, int ch2,double t1, double f1, double w1, double t2, double f2, double w2, double S, double k, double tss, double thv, int cascade){
                                                                                                                                        qodev *auxdev=(qodev *) dev;
                                                                                                                                        auxdev->add_QD(ch1,ch2,t1,f1,w1,t2,f2,w2,S,k,tss,thv,cascade); }


    void dev_add_Bell(long int dev, int ch1, int ch2,int kind, double phi, double t1, double f1, double w1, double t2, double f2, double w2){
                                                                                                                                    char ckind;
                                                                                                                                    qodev *auxdev=(qodev *) dev;
                                                                                                                                    switch (kind){
                                                                                                                                        case 0:  ckind='+'; break;
                                                                                                                                        case 1:  ckind='-'; break;
                                                                                                                                        case 2:  ckind='p'; break;
                                                                                                                                        case 3:  ckind='m'; break;
                                                                                                                                        default: ckind='+'; break;
                                                                                                                                    }
                                                                                                                                    auxdev->add_Bell(ch1,ch2,ckind,phi,t1,f1,w1,t2,f2,w2); }

    void dev_add_BellP(long int dev, int ch1, int ch2, int kind, double phi, double t1, double f1, double w1, double t2, double f2, double w2){
                                                                                                                                    char ckind;
                                                                                                                                    qodev *auxdev=(qodev *) dev;
                                                                                                                                    switch (kind){
                                                                                                                                        case 0:  ckind='+'; break;
                                                                                                                                        case 1:  ckind='-'; break;
                                                                                                                                        case 2:  ckind='p'; break;
                                                                                                                                        case 3:  ckind='m'; break;
                                                                                                                                        default: ckind='+'; break;
                                                                                                                                    }
                                                                                                                                    auxdev->add_BellP(ch1,ch2,ckind,phi,t1,f1,w1,t2,f2,w2);}
    long int dev_input(long int dev){     qodev *auxdev=(qodev *) dev; return (long int) auxdev->input();}
    long int dev_circuit(long int dev){   qodev *auxdev=(qodev *) dev; return (long int) auxdev->circuit();}
    void dev_repack(long int dev, int* ipack, int n){   qodev *auxdev=(qodev *) dev; veci tmp=to_veci(ipack,n); auxdev->repack(tmp);}
    double dev_emitted_vis(long int dev, int i, int j){ qodev *auxdev=(qodev *) dev; return auxdev->emitted_vis(i,j);}
    void dev_prnt_packets(long int dev){ qodev *auxdev=(qodev *) dev; auxdev->prnt_packets(); cout << flush;}

    // Circuit elements
    //      Basic elements
    void dev_random_circuit(long int dev){ qodev* aux=(qodev*)dev; aux->random_circuit();}
    void dev_NSX(long int  dev, int i_ch1, int i_ch2, int i_ch3){ qodev* aux=(qodev*)dev; aux->NSX(i_ch1,i_ch2,i_ch3);}
    void dev_beamsplitter(long int dev, int i_ch1, int i_ch2, double theta, double phi){qodev* aux=(qodev*)dev; aux->beamsplitter(i_ch1,i_ch2, theta, phi);}
    void dev_dielectric(long int dev, int i_ch1, int i_ch2, double ret, double imt, double rer, double imr){ cmplx t=ret+jm*imt; cmplx r=rer+jm*imr; qodev* aux=(qodev*)dev; aux->dielectric(i_ch1,i_ch2,t,r);}
    void dev_MMI2(long int dev, int i_ch1, int i_ch2){ qodev* aux=(qodev*)dev; aux->MMI2(i_ch1, i_ch2);}
    void dev_rewire(long int dev, int i_ch1, int i_ch2){ qodev* aux=(qodev*)dev; aux->rewire(i_ch1, i_ch2);}
    void dev_phase_shifter(long int dev, int i_ch, double phi){ qodev* aux=(qodev*)dev; aux->phase_shifter(i_ch,phi);}
    void dev_loss(long int dev, int i_ch, double l){ qodev* aux=(qodev*)dev; aux->loss(i_ch,l);}
    void dev_delay(long int dev, int i_ch){ qodev* aux=(qodev*)dev; aux->delay(i_ch);}

    // Polarization elements
    void dev_rotator(long int dev, int i_ch, double theta, double phi){ qodev* aux=(qodev*)dev; aux->rotator(i_ch,theta,phi);}
    void dev_pol_beamsplitter(long int dev, int i_ch1, int i_ch2, int P){qodev* aux=(qodev*)dev; aux->pol_beamsplitter(i_ch1,i_ch2, P);}
    void dev_pol_phase_shifter(long int dev, int i_ch, int P, double phi){qodev* aux=(qodev*)dev; aux->pol_phase_shifter(i_ch,P,phi);}
    void dev_pol_filter(long int dev, int i_ch, int P){qodev* aux=(qodev*)dev; aux->pol_filter(i_ch,P);}
    void dev_half(long int dev, int i_ch, double alpha){qodev* aux=(qodev*)dev; aux->half(i_ch, alpha);}
    void dev_quarter(long int dev, int i_ch, double alpha, double gamma){qodev* aux=(qodev*)dev; aux->quarter(i_ch, alpha);}

    //      Detection elements
    //      Detection elements
    void dev_detector(long int  dev, int i_ch, int cond, int pol, int mpi, int mpo, double eff, double blnk, double gamma){ qodev* aux=(qodev*)dev; aux->detector(i_ch,cond,pol,mpi,mpo,eff,blnk,gamma);}
    void dev_noise(long int dev, double stdev2){ qodev* aux=(qodev*)dev; aux->noise(stdev2);}
    long int dev_apply_condition(long int dev,long int st, bool ignore){ qodev *auxdev=(qodev *) dev; state  *auxst=(state *) st; return (long int) auxdev->apply_condition(auxst,ignore);}
    //--------------------------------------------------------------------------------------------------------------------------


    //--------------------------------------------------------------------------------------------------------------------------
    // SIMULATOR
    // Management methods
    long int sim_new_simulator(int method,int i_mem){ return (long int) new simulator(method,i_mem);}
    void sim_destroy_simulator(long int sim){simulator *aux=(simulator *)sim; delete aux; }

    // Run methods
    long int sim_run(long int sim,long int dev){ simulator *auxsim=(simulator *) sim; qodev  *auxdev=(qodev *) dev; return (long int) auxsim->run(auxdev);}
    long int sim_run_state(long int sim,long int st,long int qoc ){ simulator *auxsim=(simulator *) sim; state  *auxst=(state *) st; qocircuit *auxqoc=(qocircuit*)qoc; return (long int) auxsim->run(auxst,auxqoc);}
    // Sampling methods
    long int sim_sample(long int sim,long int dev, int N){ simulator *auxsim=(simulator *) sim; qodev  *auxdev=(qodev *) dev; return (long int) auxsim->sample(auxdev,N);}
    long int sim_sample_state(long int sim,long int st,long int qoc, int N ){ simulator *auxsim=(simulator *) sim; state  *auxst=(state *) st; qocircuit *auxqoc=(qocircuit*)qoc; return (long int) auxsim->sample(auxst,auxqoc,N);}
    //--------------------------------------------------------------------------------------------------------------------------

}
