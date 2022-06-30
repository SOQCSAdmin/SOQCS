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
    long int qoc_new_qocircuit(int i_nch, int i_nm, int i_ns, int clock, int i_R, bool loss){ return (long int)new qocircuit(i_nch,i_nm,i_ns, clock,i_R,loss);}
    void qoc_destroy_qocircuit(long int qoc){qocircuit* aux=(qocircuit*)qoc; delete aux; }
    int qoc_num_levels(long int qoc){qocircuit* aux=(qocircuit*)qoc; return aux->num_levels();}

    // Circuit elements
    //      Basic elements
    void qoc_random_circuit(long int qoc){ qocircuit* aux=(qocircuit*)qoc; aux->random_circuit();}
    void qoc_NSX(long int  qoc, int i_ch1, int i_ch2, int i_ch3){ qocircuit* aux=(qocircuit*)qoc; aux->NSX(i_ch1,i_ch2,i_ch3);}
    void qoc_beamsplitter(long int qoc, int i_ch1, int i_ch2, double theta, double phi){qocircuit* aux=(qocircuit*)qoc; aux->beamsplitter(i_ch1,i_ch2, theta, phi);}
    void qoc_dielectric(long int qoc, int i_ch1, int i_ch2, double ret, double imt, double rer, double imr){ cmplx t=ret+jm*imt; cmplx r=rer+jm*imr; qocircuit* aux=(qocircuit*)qoc; aux->dielectric(i_ch1,i_ch2,t,r);}
    void qoc_MMI2(long int qoc, int i_ch1, int i_ch2){ qocircuit* aux=(qocircuit*)qoc; aux->MMI2(i_ch1, i_ch2);}
    void qoc_phase_shifter(long int qoc, int i_ch, double ret, double imt){ cmplx t= ret+jm*imt; qocircuit* aux=(qocircuit*)qoc; aux->phase_shifter(i_ch,t);}

    //      Detection elements
    void qoc_detector(long int  qoc, int i_ch, int cond, double eff, double blnk, double gamma){ qocircuit* aux=(qocircuit*)qoc; aux->detector(i_ch,cond,eff,blnk,gamma);}
    void qoc_noise(long int qoc, double stdev2){ qocircuit* aux=(qocircuit*)qoc; aux->noise(stdev2);}

    //      Emitter and distinguishability model.
    void qoc_def_packet(long int qoc, int n, double t, double f, double w){ qocircuit* aux=(qocircuit*)qoc; aux->def_packet(n,t, f, w);}
    double qoc_emitted_vis(long int qoc, int i,int j){ qocircuit* aux=(qocircuit*)qoc; aux->emitted_vis(i, j);}
    void qoc_emitter(long int qoc, int ikind, int rand){char ckind='G';if(ikind==1) ckind='E'; qocircuit* aux=(qocircuit*)qoc; aux->emitter (ckind, rand);}
    void qoc_delay(long int qoc, int ch, double dt){ qocircuit* aux=(qocircuit*)qoc; aux->delay(ch, dt);}

    // Print functions
    void qoc_prnt(long int qoc){ qocircuit* aux=(qocircuit*)qoc; aux->prnt(); cout << flush;}
    void qoc_set_prnt_flag(qocircuit*qoc, int flag){qocircuit* aux=(qocircuit*)qoc; aux->set_prnt_flag(flag);}
    //--------------------------------------------------------------------------------------------------------------------------


    //--------------------------------------------------------------------------------------------------------------------------
    // STATE
    // Management functions
    long int st_new_state(int i_level, int i_maxket){ return (long int)new state(i_level,i_maxket);}
    void st_destroy_state(long int st){state* aux=(state*)st; delete aux; }

    // State manipulation methods.
    void st_add_term(long int st, double rampl, double iampl, int *term, int n, int m, long int qoc){ state *auxst=(state*)st;
                                                                                                      qocircuit *auxqoc=(qocircuit*)qoc;
                                                                                                      cmplx ampl=rampl+jm*iampl;
                                                                                                      mati imat=to_mati(term,n,m);
                                                                                                      auxst->add_term(ampl,imat,auxqoc);}
    long int st_post_selection(long int st, long int prj){ state* auxst=(state *)st; projector* auxprj=(projector*)prj; return (long int)auxst->post_selection(auxprj); }

    //Print methods
    void  st_prnt_state(long int st, long int qoc, bool loss,int column){  state* auxst=(state*)st; qocircuit* auxqoc=(qocircuit*)qoc; auxst->prnt_state(auxqoc,loss,column); cout << flush;}


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

    // Manipulation methods
    long int pb_calc_measure(long int pbin, long int qoc){p_bin *auxpb=(p_bin *)pbin; qocircuit *auxqoc=(qocircuit *)qoc; return (long int)auxpb->calc_measure(auxqoc);}

    //Bin consultation methods
    int pb_nbins(long int pbin){p_bin *auxpb=(p_bin*)pbin; return auxpb->nket;}
    int pb_num_levels(long int pbin){p_bin *auxpb=(p_bin*)pbin; return auxpb->nlevel;}
    char *pb_tag(long int pbin, int index){p_bin *auxpb=(p_bin*)pbin; string value =auxpb->tag(index); char* char_array=new char[value.length()+1]; strcpy(char_array, value.c_str()); return char_array;}
    double pb_prob(long int pbin,  int index){p_bin *auxpb=(p_bin*)pbin; return auxpb->prob(index);}
    double pb_prob_def(long int pbin, int *def,int n, int m, qocircuit *qoc){  p_bin *auxpb=(p_bin *) pbin;
                                                                               qocircuit *auxqoc=(qocircuit*)qoc;
                                                                               mati imat;
                                                                               imat=to_mati(def,n,m);
                                                                               auxpb->prob(imat,auxqoc);
                                                                            }
    // Print bins
    void  pb_prnt_bins(long int pbin, long int qoc, double thresh, bool loss){p_bin *auxpb=(p_bin*)pbin; qocircuit *auxqoc=(qocircuit*)qoc; auxpb->prnt_bins(auxqoc,thresh,loss); cout << flush;}
    //--------------------------------------------------------------------------------------------------------------------------


    //--------------------------------------------------------------------------------------------------------------------------
    //PHOTONS
    // Management functions
    long int ph_new_bunch(int i_level, int i_maxket){ return (long int)new ph_bunch(i_level,i_maxket);}
    void ph_destroy_bunch(long int photons){ph_bunch *aux=(ph_bunch *)photons; delete aux; }
    // Clone?

    // Manipulation methods
    void ph_add_photons(long int photons, int N, int ch, int P, double t, double f, double w, long int qoc){  ph_bunch *auxph=(ph_bunch *) photons; qocircuit *auxqoc=(qocircuit*)qoc; auxph->add_photons(N,ch,P,t,f,w,auxqoc);
                                                                                                           }
    void ph_send2circuit(long int photons, int ikind, int rand, qocircuit *qoc){char ckind='G'; if(ikind==1) ckind='E'; ph_bunch *auxph=(ph_bunch *) photons; qocircuit *auxqoc=(qocircuit*)qoc; auxph->send2circuit(ckind,rand,auxqoc);}
    void ph_weight(long int photons, double rA, double iA){ ph_bunch *auxph=(ph_bunch *) photons; cmplx A=rA+jm*iA; auxph->weight(A); }
    void ph_alternative(long int photons){ ph_bunch *auxph=(ph_bunch *) photons; auxph->alternative(); }
    long int ph_bunch_state(long int photons){ ph_bunch *auxph=(ph_bunch *) photons; return (long int) auxph->bunch_state();}

    // Print methods
    void  ph_prnt_state(long int photons, long int qoc, bool loss,int column){ ph_bunch  *auxph=(ph_bunch *) photons; qocircuit *auxqoc=(qocircuit*)qoc; auxph->prnt_state(auxqoc,loss,column); cout << flush;}
    //--------------------------------------------------------------------------------------------------------------------------


    //--------------------------------------------------------------------------------------------------------------------------
    // SIMULATOR
    // Management methods
    long int sim_new_simulator(int method,int i_mem){ return (long int) new simulator(method,i_mem);}
    void sim_destroy_simulator(long int sim){simulator *aux=(simulator *)sim; delete aux; }

    // Run methods
    long int sim_run(long int sim,long int photons,long int qoc ){ simulator *auxsim=(simulator *) sim; ph_bunch  *auxph=(ph_bunch *) photons; qocircuit *auxqoc=(qocircuit*)qoc; return (long int) auxsim->run(auxph,auxqoc);}
    long int sim_run_state(long int sim,long int st,long int qoc ){ simulator *auxsim=(simulator *) sim; state  *auxst=(state *) st; qocircuit *auxqoc=(qocircuit*)qoc; return (long int) auxsim->run(auxst,auxqoc);}
    // Sampling methods
    long int sim_sample(long int sim,long int photons,long int qoc, int N){ simulator *auxsim=(simulator *) sim; ph_bunch  *auxph=(ph_bunch *) photons; qocircuit *auxqoc=(qocircuit*)qoc; return (long int) auxsim->sample(auxph,auxqoc,N);}
    long int sim_sample_state(long int sim,long int st,long int qoc, int N ){ simulator *auxsim=(simulator *) sim; state  *auxst=(state *) st; qocircuit *auxqoc=(qocircuit*)qoc; return (long int) auxsim->sample(auxst,auxqoc,N);}
    //--------------------------------------------------------------------------------------------------------------------------
}
