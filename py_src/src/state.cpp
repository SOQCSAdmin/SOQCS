//======================================================================================================
// File state.cpp
//
// STATE DEFINITION LIBRARY.
//
// Copyright © 2022 National University of Ireland Maynooth, Maynooth University. All rights reserved.
// The contents of this file are subject to the licence terms detailed in LICENCE.TXT available in the
// root directory of this source tree. Use of the source code in this file is only permitted under the
// terms of the licence in LICENCE.TXT.
//======================================================================================================
#include "state.h"


//----------------------------------------
//
//  Creates a ket list.
//  The maximum number of kets is set by default.
//
//----------------------------------------
ket_list::ket_list(int i_level){
//  int i_level     // Number of levels to describe a ket


    create_ket_list(i_level,DEFSTATEDIM);
}


//-----------------------------------------------------
//
//  Creates a ket list specifying the maximum number of kets.
//
//-----------------------------------------------------
ket_list::ket_list(int i_level, int i_maxket){
//  int i_level     // Number of levels to describe a ket
//  int i_maxket    // Maximum number of different kets in the list


    create_ket_list(i_level,i_maxket);
}


//-----------------------------------------------------
//
//  Creates a ket list specifying the maximum number of kets
//  and the vector of equivalence between state levels and
//  circuit defined levels. Intended for internal use.
//
//-----------------------------------------------------
ket_list::ket_list(int i_level, int i_maxket, int *i_vis){
//  int  i_level     // Number of levels to describe a ket
//  int  i_maxket    // Maximum number of different kets in the list
//  int *i_vis       // Vector of equivalence between state levels and circuit defined levels
//  Auxiliary  index
    int  g;           // Aux index


    create_ket_list(i_level,i_maxket);
    for(g=0;g<i_level;g++) vis[g]=i_vis[g];
}


//----------------------------------------
//
//  Auxiliary private method to create a ket list
//  Each ket is distinguished by level occupations
//
//----------------------------------------
void ket_list::create_ket_list(int i_level, int i_maxket){
//  int i_level     // Number of levels to describe a ket
//  int i_maxket    // Maximum number of different kets in the list
//  Auxiliary index
    int i;          // Aux index


    // Init ket information
    nket=0;
    nlevel=i_level;
    maxket=i_maxket;

    // Create and int amplitude/coefficient and ket structures
    ket=new int*[maxket];
    for(i=0;i<maxket;i++){
        ket[i]=new int[nlevel]();
    }

    // Compute the trivial print visibility vector.
    vis=new int[nlevel];
    for(i=0;i<nlevel;i++){
        vis[i]=i;
    }
}


//----------------------------------------
//
//  Destroys a ket_list
//
//----------------------------------------
ket_list::~ket_list(){
//  Auxiliary index
    int i;          // Aux index


    // Liberate memory of amplitude/coefficients and ket structures
    for(i=0;i<maxket;i++){
        delete[] ket[i];
    }

    //Free memory
    delete[] ket;
    delete[] vis;
    // Clear hash table
    ketindex.clear();
}


//----------------------------------------
//
//  Copy a state
//
//----------------------------------------
ket_list *ket_list::clone(){
    // Variable
    ket_list *aux;  // Auxiliary state recipient of the copy
    // Auxiliary index
    int i;          // Aux index
    int j;          // Aux index


    aux=new ket_list(nlevel,maxket);
    aux->nket=nket;
    for(i=0;i<nket;i++){
        for(j=0;j<nlevel;j++){
            aux->ket[i][j]=ket[i][j];
        }
    }
    for(j=0;j<nlevel;j++){
        aux->vis[j]=vis[j];
    }
    aux->ketindex=ketindex;
    return aux;
}


//----------------------------------------
//
// Finds the position of a ket in the list
//
//----------------------------------------
int ket_list::find_ket(int *occ){
//  int *occ;                    // Occupation of those levels
//  Variables
    long long int value;         // Hash value in a dynamic base dictionary
    thash::const_iterator hashv; // Hash value constant iterator.


    // Update amplitude and occupation
    // Store
    value=hashval(occ,nlevel);
    hashv=ketindex.find(value);
    // If the ket does not exist yet on the output create new one.
    if(hashv==ketindex.end()){
        return -1;
    }else{
        return hashv->second;
    }
}


//----------------------------------------
//
// Finds the position of a ket in the list
// The ket is defined in human readable form
//
//----------------------------------------
int ket_list::find_ket(mati def,qocircuit *qoc){
//  mati       def;            // Definition of the ket that is searched
//  qocircuit *qoc;            // Circuit to which this list is referred
//  Variables
    ket_list  *aux;            // Bra state created from the definition def
//  Dictionary variables
    int        index;          // Index of rows


    // Create bra reference ket
    aux= new ket_list(nlevel,1,vis);
    aux->add_ket(def,qoc);

    // Search the reference
    index=find_ket(aux->ket[0]);

    // Free memory
    delete aux;

    // Return value
    return index;
}


//----------------------------------------
//
//  Adds a new ket to the list state
//
//----------------------------------------
int ket_list::add_ket(int *occ){
//  int *occ;                    // Level occupation
//  Variables
    int  index;                  // Index/List position to be returned
    long long int value;         // Hash value in a dynamic base dictionary
    thash::const_iterator hashv; // Hash value constant iterator.
//  Auxiliary index
    int  i;                      // Aux index


    // Check and warn about memory limits.
    if(nket>=maxket) cout << "Ket list add ket #1: Warning! memory limit exceeded." << endl;

    // Update amplitude and occupation
    //Store
    value=hashval(occ,nlevel);
    hashv=ketindex.find(value);

    if(hashv==ketindex.end()){
        index=nket;
        ketindex[value]=nket;
        for(i=0;i<nlevel;i++){
            ket[nket][i]=occ[i];

        }

        nket=nket+1;
    }else{
        index=hashv->second;
    }

    // Return storage index
    return index;
}


//---------------------------------------------------
//
// Adds a new ket to the list (human readable form)
//
//---------------------------------------------------
int ket_list::add_ket(hterm term, qocircuit *qoc){
//  hterm term                   // Matrix with the state specification in human readable form
//  qocircuit *qoc               // Circuit to which the the term is referred (to translate from human to numeric form)
//  Variables
    int  index;                  // Index to be returned
    int  olevel;                 // Level index
    int  ilevel;                 // Level index
    int *occ;                    // Occupation
    int *ivis;                   // Inverse visibility.
//  Auxiliary index
    int  i;                      // Aux index
    int  num;                    // Number of columns in term


    // Check and warn about memory limits.
    if(nket>=maxket) cout << "Ket list add ket #2: Warning! memory limit exceeded." << endl;

    // Initializations
    ivis=new int[qoc->nlevel];
    for(i=0;i<qoc->nlevel;i++) ivis[i]=-1;
    for(i=0;i<nlevel;i++) ivis[vis[i]]=i;
    occ=new int[nlevel]();
    for(i=0;i<nlevel;i++) occ[i]=ket[maxket-1][0];

    num=term.cols();
    // Update occupation
    for(i=0;i<num;i++){
        // The number of rows of the matrix determines its format.
        switch (term.rows()){
            case 4:     //We fully describe the state.
                olevel=qoc->i_idx[term(0,i)][term(1,i)][term(2,i)];
                ilevel=ivis[olevel];
                if(ilevel<0) cout << "Warning! Add_term: That level no longer exist in this Hilbert space" << endl;
                occ[ilevel]=term(3,i);
                break;
            case 3:     // We skip bunch (time and frequency)
                        // There are relevant in just a few problems.
                olevel=qoc->i_idx[term(0,i)][term(1,i)][0];
                ilevel=ivis[olevel];
                if(ilevel<0) cout << "Warning! Add_term: That level no longer exist in this Hilbert space" << endl;
                occ[ilevel]=term(2,i);
                break;
            case 2:     // We consider there is no polarization in this problem. (Just channel and occupation).
                olevel=qoc->i_idx[term(0,i)][0][0];
                ilevel=ivis[olevel];
                if(ilevel<0) cout << "Warning! Add_term: That level no longer exist in this Hilbert space" << endl;
                occ[ilevel]=term(1,i);
                break;
            case 1:     //  Like before we can just input occupations.
                        //  However in this case we must define all of them.
                olevel=qoc->i_idx[i][0][0];
                ilevel=ivis[olevel];
                if(ilevel<0) cout << "Warning! Add_term: That level no longer exist in this Hilbert space" << endl;
                occ[ilevel]=term(0,i);
                break;
            default:    // Something has gone wrong.
                break;
        }
    }

    // Update list
    index=add_ket(occ);

    // Free memory
    delete []occ;
    delete []ivis;

    // Return index.
    return index;
}


//-------------------------------------------------------
//
// Remove levels with wave-packet indexes greater than zero.
// Needed for some specific operations
//
//-------------------------------------------------------
ket_list *ket_list:: remove_time(qocircuit *qoc){
//  qocircuit *qoc;                 // Circuit where that is referred by this list
//  Variables
    int  newnlevel;                 // New number of levels
    int *auxket;                    // Auxiliary ket
    int *newvis;                    // New visibility vector
    int *isincluded;                // Levels included
    ket_list *auxlist;              // Auxiliary ket list (to be returned)
//  Auxiliary index
    int i;                          // Aux index
    int j;                          // Aux index
    int k;                          // Aux index


    // Init variables
    newnlevel=0;
    newvis=new int[nlevel]();
    isincluded=new int[nlevel]();
    for(i=0;i<nlevel;i++){
        if(qoc->idx[vis[i]].s==0){
                isincluded[i]=1;
                newvis[newnlevel]=vis[i];
                newnlevel=newnlevel+1;
        }
    }

    auxlist=new ket_list(newnlevel,maxket,newvis);
    delete[] newvis;

    // Compute operation
    auxket=new int[newnlevel]();
    for(i=0;i<nket;i++){
        k=0;
        for(j=0;j<nlevel;j++){
            if(isincluded[j]==1){
                auxket[k]=ket[i][j];
                k++;
            }
        }
        auxlist->add_ket(auxket);
    }

    // Free memory
    delete auxket;
    delete[] isincluded;

    // Return value
    return auxlist;
}


//-------------------------------------------------------
//
//  Clear the list (without destroying it).
//
//-------------------------------------------------------
void ket_list::clear_kets(){


    nket=0;
    ketindex.clear();
}


//--------------------------------------------
//
//  Prints a ket
//
//--------------------------------------------
void  ket_list::prnt_ket(int iket, qocircuit *qoc){
//  int iket;           // Ket number to be printed
//  qocircuit *qoc;     // Circuit to which the ket is referred. (Necessary to print in human readable form)


    prnt_ket(iket, false, qoc);
}


//--------------------------------------------
//
//  Prints a ket
//
//--------------------------------------------
void  ket_list::prnt_ket(int iket, bool loss, qocircuit *qoc){
//  int iket;           // Ket number to be printed
//  bool loss;          // Print loss channels in different color? False=No/True=Yes
//  qocircuit *qoc;     // Circuit to which the ket is referred. (Necessary to print in human readable form)
//  Variables
    int writeprev;      // Have we print a level of this ket previously? 1=Yes/0=No
    int lev;            // Level number as defined in circuit base. We do not rely in order anymore because post-selection may change it.
    int nchm;           // Half the number of channels.
//  Auxiliary index
    int k;              // Aux index


    cout  << CYAN << " | ";
    writeprev=0;
    for(k=0;k<nlevel;k++){
        lev=vis[k];
        if(qoc > nullptr){
            nchm=qoc->nch/2;
            // We print in different formats depending on the circuit flag for printing
            switch (qoc->flag_prnt){
                case 0: // Numerical form
                    if(writeprev==1) cout << ", ";
                    if((qoc->idx[lev].ch>=nchm)&&(loss)) cout  << RED;
                    cout << "( ";
                    cout << qoc->idx[lev].ch;
                    if(qoc->nm>1)    cout << ", " << qoc->idx[lev].m;
                    if(qoc->ns>1)    cout << ", " << qoc->idx[lev].s;
                    cout <<" ): " << ket[iket][k];
                    cout  << CYAN;
                    writeprev=1;
                    break;
                case 1: // Human readable form. Some numbers are changed by letters
                    if(writeprev==1) cout << ", ";
                    if((qoc->idx[lev].ch>=nchm)&&(loss)) cout  << RED;
                    cout << "( ";
                    cout << qoc->idx[lev].ch;
                    if(qoc->nm>1)    cout << ", " << pl[qoc->idx[lev].m];
                    if(qoc->ns>1)    cout << ", " << qoc->idx[lev].s;
                    cout <<" ): " << ket[iket][k];
                    cout  << CYAN;
                    writeprev=1;
                    break;
                case 2: // Human readable condensed form. A condensed form particularly apt for Baso Baset problem and similars.
                    if(ket[iket][k]>0){
                        if(writeprev==1) cout << ", ";
                        if((qoc->idx[lev].ch>=nchm)&&(loss)) cout  << RED;
                        if(ket[iket][k]>1) cout << "[" << ket[iket][k] << "]";
                        if(qoc->nm>1)    cout << pl[qoc->idx[lev].m];
                        if(qoc->ns>1)    cout << "(" << qoc->idx[lev].s << ")";
                        cout << qoc->idx[lev].ch;
                        cout  << CYAN;
                        writeprev=1;
                    }
                    break;
                case 3: // Human readable condensed form. A condensed form particularly apt for Baso Baset problem and similars.
                    if(writeprev==1) cout << ", ";
                    if((qoc->idx[lev].ch>=nchm)&&(loss)) cout  << BLUE;
                    if(ket[iket][k]>=0){
                        cout << ket[iket][k];
                    }else{
                        cout << "X";
                    }
                    cout  << CYAN;
                    writeprev=1;
                    break;

                default:
                    cout << "Prnt_state: Something went wrong.";
                    break;
            }
        }else{
            if(writeprev==1) cout << ", ";
            if(ket[iket][k]>=0){
                cout << ket[iket][k];
            }else{
                cout << "X";
            }
            writeprev=1;
        }
    }
    cout << " >"<< RESET;
}


//----------------------------------------
//
//  Creates a state
//  The maximum number of kets is set by default.
//
//----------------------------------------
state::state(int i_level):ket_list(i_level){
//  int i_level;     // Number of levels to describe the state


    ampl=new cmplx[maxket]();
}


//-----------------------------------------------------
//
//  Creates a state specifying the maximum number of kets.
//
//-----------------------------------------------------
state::state(int i_level, int i_maxket):ket_list(i_level,i_maxket){
//  int i_level     // Number of levels to describe the state
//  int i_maxket    // Maximum number of different kets in the summation


    ampl=new cmplx[maxket]();
}


//-----------------------------------------------------
//
//  Creates a state specifying the maximum number of kets
//  and the vector of equivalence between state levels and
//  circuit defined levels. Intended for internal use.
//
//-----------------------------------------------------
state::state(int i_level, int i_maxket, int *i_vis):ket_list(i_level,i_maxket, i_vis){
//  int  i_level     // Number of levels to describe the state
//  int  i_maxket    // Maximum number of different kets in the state
//  int *i_vis       // Vector of equivalence between state levels and circuit defined levels


    ampl=new cmplx[maxket]();
}


//----------------------------------------
//
//  Destroys a state
//
//----------------------------------------
state::~state(){


    // Liberate memory of amplitude
    delete[] ampl;
}


//----------------------------------------
//
//  Copy a state
//
//----------------------------------------
state *state::clone(){
    // Variable
    state *aux;  // Auxiliary state
    // Auxiliary index
    int    i;    // Aux index
    int    j;    // Aux index


    aux=new state(nlevel,maxket);
    aux->nket=nket;
    for(i=0;i<nket;i++){
        aux->ampl[i]=ampl[i];
        for(j=0;j<nlevel;j++){
            aux->ket[i][j]=ket[i][j];
        }
    }
    for(j=0;j<nlevel;j++){
        aux->vis[j]=vis[j];
    }

    aux->ketindex=ketindex;
    return aux;
}


//----------------------------------------
//
//  Clear a state
//
//----------------------------------------
void state::clear_state(){


    delete[] ampl;
    ampl=new cmplx[maxket]();
    clear_kets();
}


//----------------------------------------
//
//  Adds a new term to a state
//
//----------------------------------------
int state::add_term(cmplx i_ampl, int *occ){
//  cmplx i_ampl;       // Amplitude of the new term
//  int  *occ;          // Occupation of those levels
//  Variables
    int   index;        // Index to be returned


    index=add_ket(occ);
    ampl[index]=ampl[index]+i_ampl;

    return index;
}


//---------------------------------------------------
//
// Adds a new term to a state ( human readable form)
//
//---------------------------------------------------
int state::add_term(cmplx i_ampl,hterm term, qocircuit *qoc){
//  cmplx i_ampl;                // Amplitude of the new term
//  hterm term                   // Matrix with the state specification in human readable form
//  qocircuit *qoc               // Circuit to which the the term is referred (to translate from human to numeric form)
//  Variables
    int   index;                 // Index to be returned


    index=add_ket(term,qoc);
    ampl[index]=ampl[index]+i_ampl;

    return index;
}


//---------------------------------------------------
//
// Calculates the braket <prj|state>
//
//---------------------------------------------------
cmplx state::braket(state *bra){
//  state *bra;     // Projector with the description of the levels (and their occupations) to be post-selected
//  Variables
    cmplx  rampl;   // Amplitude to be returned
//  Auxiliary index
    int    i;       // Aux index
    int    j;       // Aux index


    rampl=0.0;
    // For each projector term
    for(i=0;i<bra->nket;i++){
            // Find the corresponding ket in the state
            j=find_ket(bra->ket[i]);
            if(j>=0){
                rampl=rampl+conj(bra->ampl[i])*ampl[j];
            }
    }

    // Return state.
    return rampl;
}


//---------------------------------------------------
//
// Normalizes the state
//
//---------------------------------------------------
void state::normalize(){
//  Variables
    double tot;     // Total/ Norm of the state
//  Auxiliary index
    int    i;       // Aux index


    tot=0.0;
    for(i=0;i<nket;i++){
        tot=tot+real(conj(ampl[i])*ampl[i]);
    }
    if(abs(tot)>xcut){
        for(i=0;i<nket;i++){
            ampl[i]=ampl[i]/sqrt(tot);
        }
    }

}

//--------------------------------------------
//
//  Post selection process of the output state
//
//--------------------------------------------
state *state::post_selection(state *prj){
//  state *prj;      // Projector with the description of the levels (and their occupations) to be post-selected
//  Variables
    int    npost;    // Number of levels in which post selection is performed
    int    selected; // Is this state selected 0=No/1=Yes
    int   *islincl;  // Is level included? 0=No/1=Yes
    int   *occ;      // Occupation
    state *nstate;   // New post-selected state
// Auxiliary index
    int    i;        // Aux index
    int    j;        // Aux index
    int    k;        // Aux index
    int    l;        // Aux index


    // Determined the number of post-selected levels and which ones
    // will be included in nstate or not.
    npost=0;
    islincl=new int[nlevel];
    for(i=0;i<prj->nlevel;i++){
        if(prj->ket[0][i]>=0){
            islincl[i]=0;
            npost++;
        }else{
            islincl[i]=1;
        }
    }

    // Create new post-selected state / Reserve memory
    nstate=new state(nlevel-npost,maxket);
    k=0;
    for(l=0;l<nlevel;l++){
        if(islincl[l]==1){
            nstate->vis[k]=vis[l];
            k++;
        }
    }

    // For each projector term
    occ=new int[nlevel-npost]();
    for(i=0;i<prj->nket;i++){
        // Post-select each ket
        for(j=0;j<nket;j++){
            // Check selection condition
            selected=1;
            k=0;
            while((k<nlevel)&&(selected==1)){
                if ((ket[j][k]!=prj->ket[i][k])&&(prj->ket[i][k]>=0))  selected=0;
                k++;
            }

            // If is selected create the list of levels and
            // occupations for those not post-selected
            if(selected==1){
                k=0;
                for(l=0;l<nlevel;l++){
                    if(islincl[l]==1){
                        occ[k]=ket[j][l];
                        k++;
                    }
                }
                nstate->add_term(ampl[j]*conj(prj->ampl[i]),occ);
            }
        }
    }


    // Free memory
    delete[] occ;
    delete[] islincl;

    // Return state.
    return nstate;
}


//-------------------------------------------------------
//
// Post select a state using a "projector" initialized to zero
// in the corresponding channels defined by ch
//
//-------------------------------------------------------
state *state:: prepare_reference(veci ch, int it, qocircuit *qoc){
//  veci       ch;                // Channel list where the post-selection is performed
//  int        it;                // Allowed "times" 0=Only 0/ 1=All
//  qocircuit *qoc;               // Circuit where the detector is placed
//  Variables
    int        nl;                // List length
    int        nch;               // Detector condition definition number of columns to read
    int        ich;               // Channel index
    int        im;                // Mode index
    int        is;                // "Time" index
    int       *ket;               // Ket
    hterm      select;            // Post-selection condition
    veci       llist;             // Level list
    projector *auxprj;            // Projector
    state     *auxstate;          // Auxiliary state
    state     *newstate;          // New post-selected state
//  Auxiliary index
    int        i;                 // Aux index
    int        l;                 // Aux index
    int        k;                 // Aux index


    // Init variables
    nch=ch.size();

    // Post select-channels
    k=0;
    if(nch>0) {
        select.setZero(4,nch*qoc->nm*qoc->ns);
        for(ich=0;ich<nch;ich++){
        for(im=0;im<qoc->nm;im++){
        for(is=0;is<qoc->ns;is++){

            select(0,k)=ch(ich);
            select(1,k)=im;
            select(2,k)=is;
            select(3,k)=0;
            k++;

        }}}

        auxprj=new projector(qoc->num_levels(),1,vis);
        auxprj->add_term(1.0,select,qoc);
        auxstate=this->post_selection(auxprj);
        delete auxprj;
    }else{
        auxstate=this->clone();
    }


    // Post-select "Times"
    if(it==0){
        nl=0;
        llist.resize(auxstate->nlevel);
        for(i=0;i<auxstate->nlevel;i++){
            l=auxstate->vis[i];
            if((qoc->idx[l].s)>0){
                llist(nl)=i;
                nl=nl+1;
            }
        }

        ket= new int[auxstate->nlevel];

        for(i=0;i<auxstate->nlevel;i++) ket[i]=-1;
        for(i=0;i<nl;i++)  ket[llist(i)]=0;

        auxprj=new projector(auxstate->nlevel,1,auxstate->vis);
        auxprj->add_term(1.0,ket);
        newstate=auxstate->post_selection(auxprj);
        delete auxprj;
        delete ket;
    }else{
        newstate=auxstate->clone();
    }

    // Free memory
    delete auxstate;

    //Return new state
    return newstate;
}


//--------------------------------------------
//
// Prints a state in numerical form
// ( various kets and amplitudes )
//
//--------------------------------------------
void  state::prnt_state(){


    prnt_state(nullptr,false,0);
}


//--------------------------------------------
//
// Prints a state in numerical form
// ( various kets and amplitudes )
//
//--------------------------------------------
void  state::prnt_state(int column){
//  int column;          // Is the state printed in line or with kets in columns? 0=No/1=Yes


    prnt_state(nullptr,false,column);
}


//--------------------------------------------
//
// Prints a state in human readable form
// ( various kets and amplitudes )
//
//--------------------------------------------
void  state::prnt_state(qocircuit *qoc, int column){
//  qocircuit *qoc;     // Circuit to which the ket is referred. (Necessary to print in human readable form)
//  int column          // Do you want the state print in line or with kets in columns? 0=No/1=Yes


    prnt_state(qoc, false, column);
}


//--------------------------------------------
//
// General method to print a state
// ( various kets and amplitudes )
//
//--------------------------------------------
void  state::prnt_state(qocircuit *qoc, bool loss,int column){
//  qocircuit *qoc;     // Circuit to which the ket is referred. (Necessary to print in human readable form)
//  bool loss;          // Print loss channels in different color? False=No/True=Yes
//  int column          // Do you want the state print in line or with kets in columns? 0=No/1=Yes


    if(column==0){
        prnt_in_rows(qoc, loss);
    }else{
        prnt_in_cols(qoc, loss);
    }

}


//--------------------------------------------
//
// General method to print a state in rows as
// a summation( various kets and amplitudes )
//
//--------------------------------------------
void  state::prnt_in_rows(qocircuit *qoc, bool loss){
//  qocircuit *qoc;     // Circuit to which the ket is referred. (Necessary to print in human readable form)
//  bool loss;          // Print loss channels in different color? False=No/True=Yes
//  Variables
    int  firstline;     // Is this the first term? 1=Yes/0=No
//  Auxiliary index
    int  i;             // Aux index


    firstline=1;
    for(i=0;i<nket;i++){
        if (abs(ampl[i])>DEFTHOLDPRNT){
            if(firstline==0) cout << " + ";

            firstline=0;

            cout << ampl[i];
            cout << " * ";
            prnt_ket(i,loss,qoc);
        }
    }
    if(firstline==1) cout << "| empty >";
    cout << endl;
}


//--------------------------------------------
//
// General method to print a state in columns as
// a list ( various kets and amplitudes )
//
//--------------------------------------------
void  state::prnt_in_cols(qocircuit *qoc, bool loss){
//  qocircuit *qoc;     // Circuit to which the ket is referred. (Necessary to print in human readable form)
//  bool loss;          // Print loss channels in different color? False=No/True=Yes
//  Variables
    int  firstline;     // Is this the first term? 1=Yes/0=No
//  Auxiliary index
    int  i;             // Aux index


    firstline=1;
    for(i=0;i<nket;i++){
        if (abs(ampl[i])>DEFTHOLDPRNT){
            prnt_ket(i,loss,qoc);


            firstline=0;
            cout << ": ";

            if(real(ampl[i])>=0){
                if(imag(ampl[i])>=0) cout << " " << fixed <<  setprecision(8) << abs(real(ampl[i])) << " + " << abs(imag(ampl[i])) << endl;
                else                 cout << " " << fixed <<  setprecision(8) << abs(real(ampl[i])) << " - " << abs(imag(ampl[i])) << endl;
            }else{
                if(imag(ampl[i])>=0) cout << "-" << fixed <<  setprecision(8) << abs(real(ampl[i])) << " + " << abs(imag(ampl[i])) << endl;
                else                 cout << "-" << fixed <<  setprecision(8) << abs(real(ampl[i])) << " - " << abs(imag(ampl[i])) << endl;
            }

        }
    }
    if(firstline==1) cout << "| empty >";
    cout << endl;
}


//---------------------------------------------------------------------
//
//  Non-ideal bell emitter with a phase e^(-i phi) in the second term.
//
//----------------------------------------------------------------------
void state::Bell(int i_ch0, int i_ch1, veci i_t,double phi,char kind, qocircuit *qoc){
//  int    i_ch0;       // Channel 0
//  int    i_ch1;       // Channel 1
//  veci   i_t;         // Vector with the wavepacket numbers to be assigned to each photon.
//  double phi;         // Non-ideal emission phase
//  char   kind;        // Kind of bell state +=Phi+/-=Phi-/p=Psi+/m=Psi-
//  qocircuit *qoc;     // Circuit to which the emitter is referred. (Necessary to interpret the defined states)
// Constants
    const  int nl=2;    // Number of level to define for the emitter
//  Variables
    double t0;          // Wavepacket 0
    double t1;          // Wavepacket 1
    double t2;          // Wavepacket 2
    double t3;          // Wavepacket 3
    complex<double> A1; // Amplitude of the first term of the bell state.
    complex<double> A2; // Amplitude of the first term of the bell state.
    mati   Bell1;       // First term of the Bell state
    mati   Bell2;       // Second term of the bell state
    state *auxemitter;  // Auxiliary emitter variable


    // Reserve memory
    Bell1.resize(DF+1,nl);
    Bell2.resize(DF+1,nl);

    // Make variables easily accessible
    t0=i_t(0);
    t1=i_t(1);
    t2=i_t(2);
    t3=i_t(3);

    // Diferent initialization depending of the kind of bell state selected
    switch(kind){
        case '+':
            Bell1 << i_ch0, i_ch1,
                         H,     H,
                         t0,     t2,
                         1,     1;
            Bell2 << i_ch0, i_ch1,
                         V,     V,
                         t1,     t3,
                         1,     1;
            A1=1.0/sqrt(2);
            A2=exp(-jm*phi)/sqrt(2.0);

        break;
        case '-':
            Bell1 << i_ch0, i_ch1,
                         H,     H,
                         t0,     t2,
                         1,     1;
            Bell2 << i_ch0, i_ch1,
                         V,     V,
                         t1,     t3,
                         1,     1;
            A1=1.0/sqrt(2);
            A2=-exp(-jm*phi)/sqrt(2.0);
        break;
        case 'p':
            Bell1 << i_ch0, i_ch1,
                         H,     V,
                         t0,     t2,
                         1,     1;
            Bell2 << i_ch0, i_ch1,
                         V,     H,
                         t1,     t3,
                         1,     1;
            A1=1.0/sqrt(2);
            A2=exp(-jm*phi)/sqrt(2.0);
        break;
        case 'm':
            Bell1 << i_ch0, i_ch1,
                         H,     V,
                         t0,    t2,
                         1,     1;
            Bell2 << i_ch0, i_ch1,
                         V,     H,
                         t1,    t3,
                         1,     1;
            A1=1.0/sqrt(2);
            A2=-exp(-jm*phi)/sqrt(2.0);
        break;
        default:
            cout << "Bell emitter: Error, undefined state";
        break;

    }


    if(nket==0){
        // If the state is not initialized just add the terms of the Bell state
        add_term( A1,Bell1,qoc);
        add_term( A2,Bell2,qoc);
    }else{
        // If the state is initialized we have to do the direct product with
        // the already existent term. We are assuming here the pre-existing
        // initialization is for different channels.
        auxemitter=new state(nlevel,maxket);
        auxemitter->add_term( A1,Bell1,qoc);
        auxemitter->add_term( A2,Bell2,qoc);
        this->dproduct(auxemitter);
        delete auxemitter;
    }
}


//---------------------------------------------------------------------
//
//  Random noise emitter
//
//----------------------------------------------------------------------
void state::Rand(int i_ch0,int i_ch1, veci i_t, qocircuit *qoc){
//  int    i_ch0;       // Channel 0
//  int    i_ch1;       // Channel 1
//  veci   i_t;         // Vector with the wavepacket numbers to be assigned to each photon.
//  qocircuit *qoc;     // Circuit to which the emitter is referred. (Necessary to interpret the defined states)
//  Constants
    const  int nl=2;    // Number of level to define for the emitter
//  Variables
    int    r;           // Random integer number between 0 and 3 (included)
    double t0;          // Wavepacket 0
    double t1;          // Wavepacket 1
    double t2;          // Wavepacket 2
    double t3;          // Wavepacket 3
    mati   Rand;        // Term of the random state
    state *auxemitter;  // Auxiliary emitter variable


    // Reserve memory
    Rand.resize(DF+1,nl);

    //Make variables easily accessible
    t0=i_t(0);
    t1=i_t(1);
    t2=i_t(2);
    t3=i_t(3);

    // Diferent initialization depending on an aleatory election
    r=(int)(floor(4.0*urand()));
    if(r==4) r=3; // Just in the improbable case we hit the limit.
    switch(r){
        case 0:
            Rand << i_ch0, i_ch1,
                         H,     H,
                         t0,     t2,
                         1,     1;

        break;
        case 1:
            Rand << i_ch0, i_ch1,
                         H,     V,
                         t0,     t3,
                         1,     1;

        break;
        case 2:
            Rand << i_ch0, i_ch1,
                         V,     H,
                         t1,     t2,
                         1,     1;

        break;
        case 3:
            Rand << i_ch0, i_ch1,
                         V,     V,
                         t1,   t3,
                         1,     1;

        break;
        default:
            cout << "Bell emitter: Error, undefined state";
        break;

    }


    if(nket==0){
        // If the state is not initialized just add the terms of the Bell state
        add_term( 1.0,Rand,qoc);
    }else{
        // If the state is initialized we have to do the direct product with
        // the already existent term. We are assuming here the pre-existing
        // initialization is for different channels.
        auxemitter=new state(nlevel,maxket);
        auxemitter->add_term(1.0,Rand,qoc);
        this->dproduct(auxemitter);
        delete auxemitter;
    }
}


//---------------------------------------------------------------------
//
//  Correlated emitter
//
//----------------------------------------------------------------------
void state::Corr(int i_ch0,int i_ch1, veci i_t, qocircuit *qoc){
//  int    i_ch0;       // Channel 0
//  int    i_ch1;       // Channel 1
//  veci   i_t;         // Vector with the wavepacket numbers to be assigned to each photon.
//  qocircuit *qoc;     // Circuit to which the emitter is referred. (Necessary to interpret the defined states)
// Constants
    const  int nl=2;    // Number of level to define for the emitter
//  Variables
    int    r;           // Random integer number between 0 and 1 (included)
    double t0;          // Wavepacket 0
    double t1;          // Wavepacket 1
    double t2;          // Wavepacket 2
    double t3;          // Wavepacket 3
    mati   Corr;        // Term of the correlated state
    state *auxemitter;  // Auxiliary emitter variable


    Corr.resize(DF+1,nl);

    t0=i_t(0);
    t1=i_t(1);
    t2=i_t(2);
    t3=i_t(3);

    // Diferent initialization depending on an aleatory election
    r=(int)(floor(2.0*urand()));
    if(r==2) r=1; // Just in the improbable case we hit the limit.
    switch(r){
        case 0:
            Corr << i_ch0, i_ch1,
                         H,     H,
                        t0,    t2,
                         1,     1;
        break;
        case 1:
            Corr << i_ch0, i_ch1,
                         V,     V,
                        t1,    t3,
                         1,     1;
        break;
        default:
            cout << "Bell emitter: Error, undefined state";
        break;

    }


    if(nket==0){
        // If the state is not initialized just add the terms of the Bell state
        add_term( 1.0,Corr,qoc);
    }else{
        // If the state is initialized we have to do the direct product with
        // the already existent term. We are assuming here the pre-existing
        // initialization is for different channels.
        auxemitter=new state(nlevel,maxket);
        auxemitter->add_term( 1.0,Corr,qoc);
        this->dproduct(auxemitter);
        delete auxemitter;
    }
}


//----------------------------------------------------------------------
//
//  Single pair photon emission in a QD
//
//----------------------------------------------------------------------
void state:: QDPair(int i_ch0,int i_ch1, veci i_t,double dt, double k, double S, double tss, double thv, qocircuit *qoc){
//  int    i_ch0;   // Channel 0
//  int    i_ch1;   // Channel 1
//  veci   i_t;     // Vector with the wavepacket numbers to be assigned to each photon.
//  double dt;      // Time between the emission of the two photons.
//  double k;       // fraction of the emitted photon pairs which originate both from a XX-X cascade in presence of background light or multiphoton emission.
//  double S;       // S/hbar FSS constant.
//  double tss;     // Time of spin-scattering tss/tx
//  double thv;     // Time of cross-dephasing thv/tx
//  qocircuit *qoc; // Circuit to which the emitter is referred. (Necessary to interpret the defined states)


    if(urand()<(k*exp(-(dt/tss)))){
        if(urand()<exp(-(dt/thv))){
            Bell(i_ch0,i_ch1,i_t,S*dt,'+',qoc);
        }else{
            Corr(i_ch0,i_ch1,i_t,qoc);
        }
    }else{
        Rand(i_ch0,i_ch1,i_t,qoc);
    }

}


//----------------------------------------------------------------------
//
//  Basso Basset Quantum Dot emitter.
//  Simplified version.
//
//----------------------------------------------------------------------
void state:: QD(mati ch, double k, double S, double tss, double thv, qocircuit *qoc){
//  mati   ch;      // Configuration matrix of the QD. First row channels. Second and third wave packet numbers of the H and V photons respectively.
//  double k;       // Fraction of the emitted photon pairs which originate both from a XX-X cascade in presence of background light or multi-photon emission.
//  double S;       // S/hbar FSS constant.
//  double tss;     // Time of spin-scattering tss/tx
//  double thv;     // Time of cross-dephasing thv/tx
//  qocircuit *qoc; // Circuit to which the emitter is referred. (Necessary to interpret the defined states)
//  Variables
    int     ch0;    // Channel 0
    int     ch1;    // Channel 1
    veci    vt;     // Wavepacket numbers vector.
    double  dt;     // Emission time t/tx.
// Note: tx is our code time unit. It is the exciton radiative lifetime.
// Auxiliary index
    int      i;     // Aux index


    vt.resize(4);
    for(i=0;i<ch.cols();i=i+2){
        // Update channel
        ch0  =ch(0,i);
        ch1  =ch(0,i+1);
        vt(0)=ch(1,i);
        vt(1)=ch(2,i);
        vt(2)=ch(1,i+1);
        vt(3)=ch(2,i+1);

        // Emit pair
        dt=-log(1-urand());
        QDPair(ch0,ch1,vt,dt,k,S,tss,thv,qoc);
    }
}


//----------------------------------------------------------------------
//
//  Direct product of states defined in non coincident channels
//  (Not to be used in general therefore private)
//
//----------------------------------------------------------------------
void state::dproduct(state *rhs){
//  state *rhs      // State in the right hand side to do the d-product
//  Variables
    int   *occ;     // Occupation
    state *aux;     // Auxiliary state
//  Auxiliary index
    int    i;       // Aux index
    int    j;       // Aux index
    int    k;       // Aux index


    // Initialize
    aux=new state(this->nlevel,this->maxket);
    aux=this->clone();
    clear_state();
    occ= new int[nlevel]();

    // Dproduct
    for(i=0;i<aux->nket;i++){
        for(j=0;j<rhs->nket;j++){
            for(k=0;k<aux->nlevel;k++){
                occ[k]=aux->ket[i][k]+rhs->ket[j][k];
            }
            add_term(aux->ampl[i]*rhs->ampl[j],occ);
        }
    }

    // Free memory
    delete[] occ;
    delete aux;
}


//-------------------------------------------------------
//
//  Creates a projector. The maximum number of kets is set by default.
//
//-------------------------------------------------------
projector::projector(int i_level):state(i_level){
//  int i_level     // Number of levels to describe the projector


    create_projector(i_level,DEFSTATEDIM);
}


//-------------------------------------------------------
//
// Creates a projector specifying the maximum number of kets.
//
//-------------------------------------------------------
projector::projector(int i_level, int i_maxket):state(i_level, i_maxket){
//  int i_level     // Number of levels to describe the projector
//  int i_maxket    // Maximum number of different kets in the summation


    create_projector(i_level,i_maxket);
}


//-----------------------------------------------------
//
//  Creates a projector specifying the maximum number of kets
//  and the vector of equivalence between state levels and
//  circuit defined levels. Intended for internal use.
//
//-----------------------------------------------------
projector::projector(int i_level, int i_maxket, int *i_vis):state(i_level, i_maxket, i_vis){
//  int i_level      // Number of levels to describe the projector
//  int  i_maxket    // Maximum number of different kets in the list
//  int *i_vis       // Vector of equivalence between state levels and circuit defined levels


    create_projector(i_level,i_maxket);
}


//-------------------------------------------------------
//
//  Auxiliary method to create a projector.
//  Same that to create a state but non-defined
//  occupation levels are not assumed to be 0 occupied.
//  They are initialized with a negative value instead.
//  This way these levels are later on ignored/not
// considered in the post-selection operation
//
//-------------------------------------------------------
void projector::create_projector(int i_level, int i_maxket){
//  int i_level     // Number of levels to describe the projector
//  int i_maxket    // Maximum number of different kets in the summation
//  Auxiliary index
    int i;          // Aux index
    int j;          // Aux index


    for(i=0;i<maxket;i++){
        for(j=0;j<nlevel;j++){
            ket[i][j]=-1;
        }
    }
}


//----------------------------------------
//
//  Creates a set of probability bins
//  The maximum number of kets is set by default.
//
//----------------------------------------
p_bin::p_bin(int i_level):ket_list(i_level){
//  int i_level     // Number of levels to describe a ket in the set


    N=0;
    p=new double[maxket]();
}


//-----------------------------------------------------
//
//  Creates a set of probability bins specifying the maximum number of bins.
//
//-----------------------------------------------------
p_bin::p_bin(int i_level, int i_maxket):ket_list(i_level,i_maxket){
//  int i_level     // Number of levels to describe a ket in the set
//  int i_maxket    // Maximum number of different kets in the set


    N=0;
    p=new double[maxket]();
}


//-----------------------------------------------------
//
//  Creates a set of probability bins specifying the maximum
//  number of bins  and the vector of equivalence between bin
//  levels and circuit defined levels. Intended for internal use.
//
//-----------------------------------------------------
p_bin::p_bin(int i_level, int i_maxket, int *i_vis):ket_list(i_level,i_maxket, i_vis){
//  int i_level     // Number of levels to describe a ket in the set
//  int i_maxket    // Maximum number of different kets in the set
//  int *i_vis      // Vector of equivalence between state levels and circuit defined levels


    N=0;
    p=new double[maxket]();;
}


//----------------------------------------
//
//  Destroys a state
//
//----------------------------------------
p_bin::~p_bin(){


    // Liberate memory of amplitude
    delete[] p;
}


//----------------------------------------
//
//  Copy a state
//
//----------------------------------------
p_bin *p_bin::clone(){
    // Variable
    p_bin *aux;     // Auxiliary state
    // Auxiliary index
    int    i;       // Aux index
    int    j;       // Aux index


    aux=new p_bin(nlevel,maxket);
    aux->nket=nket;
    for(i=0;i<nket;i++){
        aux->p[i]=p[i];
        for(j=0;j<nlevel;j++){
            aux->ket[i][j]=ket[i][j];
        }
    }
    for(j=0;j<nlevel;j++){
        aux->vis[j]=vis[j];
    }

    aux->ketindex=ketindex;
    aux->N=N;
    return aux;
}


//----------------------------------------
//
//  Counts a new sample
//
//----------------------------------------
int p_bin::add_count(int *occ){
//  int *occ;           // Occupation of those levels
//  Variables
    int index;          // Index to be returned


    index=add_ket(occ);
    p[index]=p[index]+1.0;
    N=N+1;

    return index;
}


//----------------------------------------
//
//  Adds/sums the statistic from another
//  compatible bin set.
//
//----------------------------------------
void p_bin::add_bin(p_bin *input){
//  p_bin *input;     // Input probability bin set
//  Variables
    int    index;     // Index of each entry in this set
//  Auxiliary index
    int    i;         // Aux index


    for(i=0;i<input->nket;i++){
        index=add_ket(input->ket[i]);
        p[index]=p[index]+input->p[i];
    }

    N=N+input->N;
}


//----------------------------------------
//
// Adds/ sums the statistics form a quantum
// state computing the probability outcomes
// from it
//
//----------------------------------------
void p_bin::add_state(state *input){
//  state *input;     // Input state
//  Variables
    int    index;     // Index of each entry in this set
//  Auxiliary index
    int    i;         // Aux index


    for(i=0;i<input->nket;i++){
        index=add_ket(input->ket[i]);
        p[index]=p[index]+(double) abs(conj(input->ampl[i])*input->ampl[i]);
    }

    N=N+1;
}


//----------------------------------------
//
// Returns the occupation of the indexed bin
// in string format.
//
//----------------------------------------
string p_bin::tag(int index){
//  int index;             // Index of the bin
//  Variables.
    int len;               // Length of the auxiliary string
    int diff;              // Difference with the number of levels
    long long int value;   // Hash value of the ket in that idex.
    string aux;            // Auxiliary string. Conversion of value
    string temp;           // Temporary string=aux+0's at the left to match nlevel length
//  Auxiliary index
    int i;                 // Aux index


    // Perform conversion
    value=hashval(ket[index],nlevel,9);
    aux=to_string(value);
    // Match the nlevel lenfth with 0's at the left
    len=aux.length();
    diff=nlevel-len;
    for(i=0;i<diff;i++) temp.append("0",1);
    temp.append(aux);
    // Return result
    return temp;
}


//----------------------------------------
//
// Probability value of the indexed bin
//
//----------------------------------------
double p_bin::prob(int index){
//  int index;            // Index of the bin


    return p[index]/(double)N;
}


//----------------------------------------
//
// Probability value a bin from it definition
//
//----------------------------------------
double p_bin::prob(mati def,qocircuit *qoc){
//  mati       def;          // Definition of a probability bin in human readable form
//  qocircuit *qoc;          // Circuit where the detector is placed
//  Variables
    int index;               // Index in this set
    ket_list  *bra;          // Bra state created from the definition def


    // Create bra state
    bra= new ket_list(nlevel,1,vis);
    bra->add_ket(def,qoc);
    index=find_ket(bra->ket[0]);
    delete bra;

    // Return probability value if present.
    if(index>=0) return prob(index);
    return 0.0;
}


//--------------------------------------------
//
// Prints the bin list
// ( various entries and probabilities )
//
//--------------------------------------------
void  p_bin::prnt_bins(){


    prnt_bins(nullptr,xcut,false);
}

//--------------------------------------------
//
// Prints the bin list
// Values below thresh are not printed
//
//--------------------------------------------
void  p_bin::prnt_bins(qocircuit *qoc, double thresh){
//  qocircuit *qoc;     // Circuit to which the ket is referred. (Necessary to print in human readable form)
//  double thresh;      // Probabilities below this value are not printed


    prnt_bins(qoc,thresh,false);
}


//--------------------------------------------
//
// Prints the bin list
// ( in human readable form )
//
//--------------------------------------------
void p_bin::prnt_bins(qocircuit *qoc, double thresh, bool loss){
//  qocircuit *qoc;     // Circuit to which the ket is referred. (Necessary to print in human readable form)
//  double thresh;      // Probabilities below this value are not printed
//  bool   loss;        // Print loss channels in different color? False=No/True=Yes
//  Variables
    int    firstline;   // Is this the first term? 1=Yes/0=No
    double prob;        // Probability
//  Auxiliary index
    int  i;             // Aux index


    firstline=1;
    for(i=0;i<nket;i++){
        prob= p[i]/(double) N;
        if (prob>thresh){
            if(firstline==0){
                cout << endl;
            }
            firstline=0;

            cout <<right << setw(2) << i << " : ";
            prnt_ket(i,loss,qoc);
            cout << ": ";
            cout << left << prob;


        }
    }
    if(firstline==1) cout << "| empty >";
    cout << endl;
}


//--------------------------------------------
//
//  Post selection process of the output state
//
//--------------------------------------------
p_bin *p_bin::post_selection(state *prj){
//  state *prj;       // Projector with the description of the levels (and their occupations) to be post-selected
//  Variables
    int    index;     // Index of the sampled state
    int    npost;     // Number of levels in which post selection is performed
    int    selected;  // Is this state selected 0=No/1=Yes
    int   *islincl;   // Is level included? 0=No/1=Yes
    int   *occ;       // Occupation
    p_bin *nbin;      // New post-selected set of probabilities
// Auxiliary index
    int    i;         // Aux index
    int    j;         // Aux index
    int    k;         // Aux index
    int    l;         // Aux index


    // Determined the number of post-selected levels and which ones
    // will be included in nstate or not.
    npost=0;
    islincl=new int[nlevel];
    for(i=0;i<prj->nlevel;i++){
        if(prj->ket[0][i]>=0){
            islincl[i]=0;
            npost++;
        }else{
            islincl[i]=1;
        }
    }

    // Create new post-selected state / Reserve memory
    nbin=new p_bin(nlevel-npost,maxket);
    k=0;
    for(l=0;l<nlevel;l++){
        if(islincl[l]==1){
            nbin->vis[k]=vis[l];
            k++;
        }
    }

    // For each projector term
    occ=new int[nlevel-npost]();
    for(i=0;i<prj->nket;i++){
        // Post-select each ket
        for(j=0;j<nket;j++){
            // Check selection condition
            selected=1;
            k=0;
            while((k<nlevel)&&(selected==1)){
                if ((ket[j][k]!=prj->ket[i][k])&&(prj->ket[i][k]>=0))  selected=0;
                k++;
            }

            // If is selected create the list of levels and
            // occupations for those not post-selected
            if(selected==1){
                k=0;
                for(l=0;l<nlevel;l++){
                    if(islincl[l]==1){
                        //nstate->ket[nstate->nket][k]=ket[j][l];
                        occ[k]=ket[j][l];
                        k++;
                    }
                }
                index=nbin->add_ket(occ);
                nbin->p[index]=nbin->p[index]+p[j]*pow(abs(prj->ampl[i]),2);
            }
        }
    }
    nbin->N=N;

    // Free memory
    delete[] occ;
    delete[] islincl;

    // Return bin
    return nbin;
}


//--------------------------------------------
//
//  Calculate a measurement from the stored statistics
//  using the corresponding circuit detector definitions
//
//--------------------------------------------
p_bin *p_bin::calc_measure(qocircuit *qoc){
//  qocircuit *qoc;         // Circuit with the detector definitions to calculate the measurement.
//  Variables
    int    S;               // Blinking and dark counts calculation number of iterations
    double stdev;           // Gaussian white noise standard deviation
    p_bin *dark;            // Output after adding dark counts
    p_bin *blinked;         // Output after blinking effect considerations
    p_bin *measured;        // Output after considering detector conditions and removing loss channels (if needed)
    p_bin *counted;         // Output after removing temporal degrees of freedom (if needed)
    p_bin *noisy;           // Output after adding some gaussian white noise
    p_bin *aux;             // Auxiliary probability bin


    // Itialize variables
    S=qoc->R;
    stdev=sqrt(qoc->dev);

    // Compute Dark counts
    if(qoc->timed==0) dark=this->dark_counts(S,qoc);
    else dark=this->clone();

    // Compute Blinking
    blinked=dark->blink(S,qoc);

    // Measure
    if(qoc->losses==1) measured=blinked->compute_loss(qoc->ncond,qoc->det_def,qoc);
    else measured=blinked->compute_cond(qoc->ncond,qoc->det_def,qoc);

    // Remove time
    if(qoc->ns>1){
       switch(qoc->timed){
            case 0:  // Counter
                aux=measured->perform_count(qoc);
                counted=aux->remove_time(qoc);
                delete aux;
                break;
            case 1:  // Clock: Partial trace
                counted=measured->remove_freq(qoc);
                break;
            case 3:  // Clock: Manual mode
                counted=measured->remove_freq(qoc);
                break;
            case 2:  // Clock+Spectrum: Full
                counted=measured->clone();
                break;
            default: // Default: Error
                cout << "P_Bin calc_measure: Time has a non-valid value" << endl;
                counted=this->clone();
                break;

       }
    }else{
        counted=measured->clone();

    }


    // Add noyse
    if (stdev>xcut) noisy=counted->white_noise(stdev);
    else noisy=counted->clone();

    // Free memory
    delete dark;
    delete blinked;
    delete measured;
    delete counted;

    // Return final detection
    return noisy;
}


//--------------------------------------------
//
// Adds Gaussian white noise to the stored statistics
//
//--------------------------------------------
p_bin *p_bin::white_noise(double stdev){
//  double stdev;       // Standard deviation of the noise
//  Variables
    double value;       // Sampled value of the normal distribution
    p_bin *newpbin;     // New output probability bin
//  Auxiliary index
    int    i;           // Aux index


    newpbin=this->clone();
    for(i=0;i<newpbin->nket;i++){
        value=grand(0.0,stdev);
        newpbin->p[i]=newpbin->p[i]+newpbin->N*value;
        if (newpbin->p[i]<0) newpbin->p[i]=0.0;
    }

    // Return new bin
    return newpbin;
}


//--------------------------------------------
//
// Computes the effect of the detector blinking
//
//--------------------------------------------
p_bin *p_bin::blink(int S, qocircuit* qoc){
//  int S;              // Number of iterations to perform the computation
//  qocircuit *qoc;     // Circuit where the detectors are defined
//  Variables
    int    ch;          // Channel number
    int    m;           // Mode number
    int    s;           // Wavepacket number
    int    blnk;        // Blink channel 0=No/1=yes
    int    index;       // Index of a ket in this set of bins
    int   *occ;         // Level occupation
    double prob;        // Blinking probability
    p_bin *newpbin;     // New output probability bin
//  Auxiliary index
    int    i;           // Aux index
    int    j;           // Aux index
    int    k;           // Aux index
    int    l;           // Aux index


    // If the number of iterations is zero just return the input
    if(S==0){
        newpbin=this->clone();
        return newpbin;
    }

    // Calculate blinking
    newpbin=new p_bin(nlevel,maxket,vis);
    for(j=0;j<S;j++){
        for(i=0;i<nket;i++){
            occ=new int[nlevel]();
            for(s=0;s<qoc->ns;s++){
            for(ch=0;ch<qoc->ndetc;ch++){
                // For each channel calculate the probability of blinking at some determined time
                prob=urand();
                if(prob>=qoc->det_par(0,ch)) blnk=0;
                else blnk=1;

                // Store the level occupation depending if the channel is blinking at that moment
                for(m=0;m<qoc->nm;m++){
                    l=qoc->i_idx[ch][m][s];
                    k=0;
                    while(vis[k]!=l) k=k+1;

                    if(blnk==0) occ[k]=occ[k]+ket[i][k];
                    else occ[k]=0;
                }
            }}

            //Loss channel do not blink. Their occupations are just copied
            for(s=0;s<qoc->ns;s++){
            for(ch=qoc->ndetc;ch<qoc->nch;ch++){
                for(m=0;m<qoc->nm;m++){
                    l=qoc->i_idx[ch][m][s];
                    k=0;
                    while(vis[k]!=l) k=k+1;

                    occ[k]=occ[k]+ket[i][k];
                }
            }}

            // Store level occupation
            index=newpbin->add_ket(occ);
            newpbin->p[index]=newpbin->p[index]+p[i];
            delete occ;
        }
    }
    newpbin->N=S*N;

    // Return new bin
    return newpbin;
}


//--------------------------------------------
//
// Computes the effect of the addition of dark counts
//
//--------------------------------------------
p_bin *p_bin::dark_counts(int S, qocircuit* qoc){
//  int S;              // Number of iterations to perform the computation
//  qocircuit *qoc;     // Circuit where the detectors are defined
//  Variables
    int    ch;          // Channel number
    int    m;           // Mode number
    int    number;      // Number of photons created by dark counts
    int    total;       // Total number of photons create by dark counts
    int    nonzero;     // Number of channel with non-zero additions due to dark counts
    int    index;       // Index of a ket in this set of bins
    int   *occ;         // Level occupation
    p_bin *newpbin;     // New output probability bin
//  Auxiliary index
    int    i;           // Aux index
    int    k;           // Aux index
    int    l;           // Aux index


    // If the number of iterations is zero just return the input
    newpbin=this->clone();
    if(S==0) return newpbin;

    // Calculate Dark counts
    nonzero=0;
    for(i=0;i<newpbin->nket;i++) newpbin->p[i]=S*newpbin->p[i];
    for(i=0;i<S;i++){
        occ=new int[nlevel]();
        total=0;
        for(ch=0;ch<qoc->ndetc;ch++){
        for(m=0;m<qoc->nm;m++){
                // For each channels the number of extra dark counts photons is obtained
                number=prand(qoc->det_par(1,ch));
                l=qoc->i_idx[ch][m][0];
                k=0;
                while(vis[k]!=l) k=k+1;

                occ[k]=number;
                total=total+number;
        }}
        // If we add some new photon the new ket is stored
        if(total>0){
            index=newpbin->add_ket(occ);
            newpbin->p[index]=newpbin->p[index]+(double)N;
            nonzero=nonzero+1;
        }
        delete occ;
    }
    newpbin->N=newpbin->N*S;

    // Return new bin
    return newpbin;
}


//----------------------------------------
//
// Computes the effect of conditional detection
//
//----------------------------------------
p_bin *p_bin:: compute_cond(int ndec, mati def,qocircuit *qoc){
//  int ndec;                       // Number of detector. ndec<nch because loss channels have no detector
//  mati       def;                 // Definition of the conditions of detection
//  qocircuit *qoc;                 // Circuit where the detector is placed
//  Variables
    int isnotfinished1;             // Is finished polarization loop? 1=Yes/0=No
    int isnotfinished2;             // Is finished time loop? 1=Yes/0=No
    hterm select;                   // Post-selection condition
    p_bin *post_selected;           // New post-selected state for a single projector
    p_bin *conditioned;             // New post-selected state
    projector *prj;                 // Projector
    int *tim;                       // "Time" configuration vector
    int *pol;                       // Polarization configuration vector
    int key[3];                     // Key of the projector definition
    int *keyprj;                    // Key of the projector
    long int selvalue;              // Projector definition index
    long int prjvalue;              // Projector index
    thash  selhash;                 // Definitions hash table
    thash  prjhash;                 // Projectors hash table.
    thash::const_iterator vselhash; // Row hash value.
    thash::const_iterator vprjhash; // Row hash value.
    veci ch;                        // List of channels by photons
    int nph;                        // Number of photons
    int maxch;                      // Larger channel number
    int selbase;                    // Base of the projector definition hash table
    int nempty;                     // Number of channels with zero fotons
    int nentry;                     // Number of entries in the projector definition
    int nprj;                       // Number of projectors
    int iph;                        // Index of photons
    int ich;                        // Channel index
    int im;                         // Mode index
    int is;                         // "Time" index
//  Auxiliary index
    int i;                          // Aux index.
    int j;                          // Aux index
    int k;                          // Aux index
    int l;                          // Aux index
    int first;


    if(ndec>0){
    conditioned=new p_bin(nlevel,maxket,vis);
    first=1;

    // Calculate:
    // How many photons.
    // How many empty channels
    // Maximum number of photons in channel
    nph=0;
    nempty=0;
    for(ich=0;ich<ndec;ich++){
        nph=nph+def(1,ich);
        maxch=max(maxch,def(0,ich));
        if(def(1,ich)==0) nempty=nempty+1;

    }
    selbase=max(maxch,max(qoc->nm,qoc->ns));

    // Create a list for every photon in which channel they are.
    // Empty channels that are post-selected to zero are appended at the endl.
    ch.resize(nph+nempty);
    k=0;
    l=0;
    for(ich=0;ich<ndec;ich++){
        for(iph=0;iph<def(1,ich);iph++){
            ch[k]=def(0,ich);
            k++;
        }

        if(def(1,ich)==0){
            ch(nph+l)=def(0,ich);
            l++;
        }
    }


    // Initialize projector aux variables.
    nprj=0;
    keyprj=new int[ndec*qoc->nm*qoc->ns]();

    // From the detector definition of conditional detection create
    // all the projections necessary considering all the possible
    // polarizations and times.
    // LOOP Polarization¡
    isnotfinished1=1;
    if(nph>0) pol=new int[nph]();
    else pol=new int[1]();

    while(isnotfinished1==1){
        // LOOP "Time"
        isnotfinished2=1;
        if(nph>0) tim=new int[nph]();
        else tim=new int[1]();

        while(isnotfinished2==1){
            // Create Projector definition
            nentry=0;
            selhash.clear();
            select.setZero(4,ndec*qoc->nm*qoc->ns);
            for(iph=0;iph<(nph+nempty);iph++){
            for(im=0;im<qoc->nm;im++){
            for(is=0;is<qoc->ns;is++){
                key[0]=ch(iph);
                key[1]=im;
                key[2]=is;


                // We have to check if for those photons the entry in select
                // already exists. If exist we update otherwise it has to be
                // updated. Otherwise we create it.
                selvalue=hashval(key,3,selbase);
                vselhash=selhash.find(selvalue);
                if(vselhash==selhash.end()){
                    k=nentry;
                    selhash[selvalue]=nentry;
                    nentry=nentry+1;
                }else{
                    k=vselhash->second;
                }

                select(0,k)=ch(iph);
                select(1,k)=im;
                select(2,k)=is;

                if(iph<nph){
                    if((im==pol[iph])&&(is==tim[iph])) select(3,k)=select(3,k)+1;
                }
            }}}

            // Create projector.
            // Check that we have not generated that projector before.
            // Photons are indistinguishable therefore in the way this loop
            // is built repetitions may appear. If it has not been generated
            // before the projector is created
            for(i=0;i<(ndec*qoc->nm*qoc->ns);i++){
                keyprj[i]=select(3,i);
            }
            prjvalue=hashval(keyprj,ndec*qoc->nm*qoc->ns);
            vprjhash=prjhash.find(prjvalue);
            // If the projector has not been created before
            if(vprjhash==prjhash.end()){
                // We updated the entry in the corresponding hash table
                prjhash[prjvalue]=nprj;
                nprj=nprj+1;

                // We create the projector
                prj=new projector(qoc->num_levels(),1,vis); //<<<<Check this
                prj->add_term(1.0,select,qoc);

                // Apply the projector to the input state and store the result
                post_selected=this->post_selection(prj);
                if(first==1){
                    delete conditioned;
                    conditioned=new p_bin(post_selected->nlevel,maxket,post_selected->vis);
                    first=0;
                }
                conditioned->add_bin(post_selected);

                // Free memory
                delete post_selected;
                delete prj;
            }

            // ADVANCE "Time" in LOOP
            if(nph>0){
                tim[0]=tim[0]+1;
                for(j=0;j<nph;j++){
                    if(tim[j]>=qoc->ns){
                        tim[j]=0;
                        if((j+1)<nph) tim[j+1]=tim[j+1]+1;
                        else isnotfinished2=0;
                    }
                }
            }else{
                isnotfinished2=0;
            }
        }

        // ADVANCE Polarization in LOOP
        if(nph>0){
            pol[0]=pol[0]+1;
            for(j=0;j<nph;j++){
                if(pol[j]>=qoc->nm){
                    pol[j]=0;
                    if((j+1)<nph) pol[j+1]=pol[j+1]+1;
                    else isnotfinished1=0;
                }
            }
        }else{
            isnotfinished1=0;
        }
        // At the end of the loop free memory
        delete[] tim;
    }
    // At the end of the loop free memory
    delete[] pol;

    // Free memory;
    delete[] keyprj;
    selhash.clear();
    prjhash.clear();
    }else{
        conditioned=this->clone();
    }

    conditioned->N=N;

    // Return new bin
    return conditioned;
}


//----------------------------------------
//
// Computes the effect of detection considering losses
//
//----------------------------------------
p_bin *p_bin:: compute_loss(int ndec, mati def, qocircuit *qoc){
//  int ndec;                       // Number of detector. ndec<nch because loss channels have no detector
//  mati       def;                 // Definition of the conditions of detection
//  qocircuit *qoc;                 // Circuit where the detector is placed
//  Variables
    int nloss;                      // Number of photons lost;
    int fchloss;                    // First loss channel
    int nchloss;                    // Number of loss channels
    int nsel;                       // Number of post-selection conditions
    int first;                      // First bin calculated 1=Yes/0=no
    int isnotfinished;              // Is not finished flag
    p_bin *conditioned;             // New post-selected state for a single projector
    p_bin *lossed;                  // New post-selected state
    veci chloss;                    // In which channel is each lost photon
    veci occ;                       // Occupation of the channels with losses
    mati select;                    // Post-selection condition of the losses
    mati selectfull;                // Full post-selection condition
//  Hash table
    int *keysel;                    // Key of the projector
    long int selvalue;              // Projector definition index
    thash  selhash;                 // Definitions hash table
    thash::const_iterator vselhash; // Row hash value.
//  Auxiliary index
    int i;                          // Aux index
    int j;                          // Aux index
    int k;                          // Aux index



    lossed=new p_bin(nlevel,maxket,vis);
    first=1;

    // Init variables
    nchloss=qoc->nch/2;
    fchloss=qoc->nch/2;

    // Reserve memory
    keysel=new int(nchloss);
    // For any number of photon loss
    for(nloss=0;nloss<=base;nloss++){
        selhash.clear();
        nsel=0;
        // We initialize in which channel is each photon
        if(nloss>0) chloss.setZero(nloss);

        // We post_select all the configurations of nloss photons
        isnotfinished=1;
        while(isnotfinished){
            // First calculate the occupations
            occ.setZero(nchloss);
            for(i=0;i<nloss;i++) occ(chloss(i))=occ(chloss(i))+1;

            select.resize(2,nchloss);
            for(k=0;k<nchloss;k++){
                select(0,k)=fchloss+k;
                select(1,k)=occ(k);
                keysel[k]=occ(k);
            }


            // Do post-selection but check if repeated
            selvalue=hashval(keysel,nchloss,nloss);
            vselhash=selhash.find(selvalue);
            // If the projector has not been created before
            if(vselhash==selhash.end()){
                // We updated the entry in the corresponding hash table
                selhash[selvalue]=nsel;
                nsel=nsel+1;
                selectfull.resize(2,select.cols()+ndec);
                k=0;
                for(i=0;i<ndec;i++){
                    selectfull(0,k)=def(0,i);
                    selectfull(1,k)=def(1,i);
                    k=k+1;
                }
                for(i=0;i<select.cols();i++){
                    selectfull(0,k)=select(0,i);
                    selectfull(1,k)=select(1,i);
                    k=k+1;
                }

                conditioned=this->compute_cond(selectfull.cols(),selectfull,qoc);
                if(first==1){
                    delete lossed;
                    lossed=new p_bin(conditioned->nlevel,maxket,conditioned->vis);
                    first=0;
                }
                lossed->add_bin(conditioned);

                delete conditioned;
            }


            // Update in which channel is each photon to sweep them all.
            if(nloss>0){
                chloss(0)=chloss(0)+1;
                for(j=0;j<nloss;j++){
                    if(chloss(j)>=nchloss){
                        chloss(j)=0;
                        if((j+1)<nloss) chloss(j+1)=chloss(j+1)+1;
                        else isnotfinished=0;
                    }
                }
            }else{
                isnotfinished=0;
            }

        }
    }

    // Properly update measured state number
    lossed->N=N;

    //Free memory
    delete[] keysel;

    // Return new bin
    return lossed;

}


//----------------------------------------
//
// Remove frequencies and leaves only times
//
//----------------------------------------
p_bin *p_bin::remove_freq(qocircuit* qoc){
//  qocircuit *qoc;     // Circuit where the detectors are defined
//  Variables
    int    ch;          // Channel number
    int    m;           // Mode number
    int    s;           // Packet number
    int    index;       // Index of a ket in this set of bins
    int   *occ;         // Level occupation
    p_bin *newpbin;     // New output probability bin
//  Auxiliary index
    int    i;           // Aux index
    int    j;           // Aux index
    int    k;           // Aux index
    int    l;           // Aux index


    // If for some reason there are not packets defined
    // then the input is returned
    if(qoc->emitted->pack_def.cols()==0){
        newpbin=this->clone();
        return newpbin;
    }

    // Calculate counts in a channel-polarization
    newpbin=new p_bin(nlevel,nket,vis);
    for(i=0;i<nket;i++){
        occ=new int[nlevel]();
        for(j=0;j<nlevel;j++){
            // We put all the counts in
            // the same time/wavepacket
            ch=qoc->idx[vis[j]].ch;
            m=qoc->idx[vis[j]].m;
            s=qoc->idx[vis[j]].s;

            if(ket[i][j]>0){ // We prevent to check levels that we don't have packet defined when we reserve more packets than we use
                // Redefine the index as time instead of packets
                l=qoc->i_idx[ch][m][qoc->emitted->pack_def(0,s)];

                // Store
                k=0;
                while(vis[k]!=l) k=k+1;
                occ[k]=occ[k]+ket[i][j];
            }
        }
        index=newpbin->add_ket(occ);
        newpbin->p[index]=newpbin->p[index]+p[i];
        delete occ;
    }
    newpbin->N=N;

    // Return new bin.
    return newpbin;
}


//----------------------------------------
//
// Performs the count of photons in a channel
// independently of time and frequency
//
//----------------------------------------
p_bin *p_bin::perform_count(qocircuit* qoc){
//  qocircuit *qoc;     // Circuit where the detectors are defined
//  Variables
    int    ch;          // Channel number
    int    m;           // Mode number
    int    index;       // Index of a ket in this set of bins
    int   *occ;         // Level occupation
    p_bin *newpbin;     // New output probability bin
//  Auxiliary index
    int    i;           // Aux index
    int    j;           // Aux index
    int    k;           // Aux index
    int    l;           // Aux index


    // Calculate counts in a channel-polarization
    newpbin=new p_bin(nlevel,nket,vis);
    for(i=0;i<nket;i++){
        occ=new int[nlevel]();
        for(j=0;j<nlevel;j++){
            // We put all the counts in
            // the same time/wavepacket
            ch=qoc->idx[vis[j]].ch;
            m=qoc->idx[vis[j]].m;
            l=qoc->i_idx[ch][m][0];

            // Store
            k=0;
            while(vis[k]!=l) k=k+1;
            occ[k]=occ[k]+ket[i][j];
        }
        index=newpbin->add_ket(occ);
        newpbin->p[index]=newpbin->p[index]+p[i];
        delete occ;
    }
    newpbin->N=N;

    // Return new bin.
    return newpbin;
}


//----------------------------------------
//
// Remove all the packet degrees of freedom
// except the 0.
//
//----------------------------------------
p_bin *p_bin:: remove_time(qocircuit *qoc){
//  qocircuit *qoc;              // Circuit where the detector is placed
//  Variables
    int    index;                // Index of a ket in this set of bins
    int    newnlevel;            // New number of levels
    int   *auxket;               // Auxiliary ket
    int   *newvis;               // New visibility vector
    int   *isincluded;           // Levels included
    p_bin *auxlist;              // Auxiliary ket list (to be returned)
//  Auxiliary index
    int    i;                    // Aux index
    int    j;                    // Aux index
    int    k;                    // Aux index


    // Init variables
    newnlevel=0;
    newvis=new int[nlevel]();
    isincluded=new int[nlevel]();
    for(i=0;i<nlevel;i++){
        if(qoc->idx[vis[i]].s==0){
                isincluded[i]=1;
                newvis[newnlevel]=vis[i];
                newnlevel=newnlevel+1;
        }
    }
    auxlist=new p_bin(newnlevel,maxket,newvis);
    delete[] newvis;

    // Perform operations
    auxket=new int[newnlevel]();
    for(i=0;i<nket;i++){
        k=0;
        for(j=0;j<nlevel;j++){
            if(isincluded[j]==1){
                auxket[k]=ket[i][j];
                k++;
            }
        }
        index=auxlist->add_ket(auxket);
        auxlist->p[index]=auxlist->p[index]+p[i];
    }
    auxlist->N=N;

    // Free memory
    delete auxket;
    delete[] isincluded;

    // Return new bin
    return auxlist;
}


//----------------------------------------
//
//  Creates a photon bunch
//  The maximum number of kets is set by default.
//
//----------------------------------------
ph_bunch::ph_bunch(int i_level){
//  int i_level;     // Number of levels to describe the state


    create_ph_bunch(i_level,DEFSTATEDIM);
}


//-----------------------------------------------------
//
//  Creates a photon bunch specifying the maximum number of kets.
//
//-----------------------------------------------------
ph_bunch::ph_bunch(int i_level, int i_maxket){
//  int i_level;     // Number of levels to describe the state
//  int i_maxket;    // Maximum number of bunches in the list


    create_ph_bunch(i_level,i_maxket);
}


//----------------------------------------
//
//  Auxiliary private method to create a photon bunch
//
//----------------------------------------
void ph_bunch::create_ph_bunch(int i_level, int i_maxket){
//  int  i_level;     // Number of levels to describe the state
//  int  i_maxket;    // Maximum number of bunches in the list
//  Variables
    int *occ;         // Level occupation

    npack=0;
    pack_list.resize(4,base);
    bunch=new state(i_level,i_maxket);
    occ=new int[i_level]();
    bunch->add_ket(occ);

    delete occ;
}


//----------------------------------------
//
//  Destroys a photon bunch object
//
//----------------------------------------
ph_bunch::~ph_bunch(){


    delete bunch;

}


//----------------------------------------
//
//  Clears a photon bunch
//
//----------------------------------------
void ph_bunch::clear(){
//  int  i_level;     // Number of levels to describe the state
//  int  i_maxket;    // Maximum number of bunches in the list
//  Variables
    int   *occ;       // Level occupation
    state *newbunch;  // New state of the photons

    npack=0;
    newbunch=new state(bunch->nlevel,bunch->maxket);
    delete bunch;
    bunch=newbunch;
    occ=new int[bunch->nlevel]();
    bunch->add_ket(occ);

    delete occ;
}


//----------------------------------------
//
//  Add photons to a photon bunch
//
//----------------------------------------
void ph_bunch::add_photons(int N, int ch, qocircuit *qoc){
//  int N;      // Number of photons to be added
//  int ch;     // Chanel where the photons are added
//  qocircuit* qoc; // Circuit to which we will emit the photons


    add_photons(N,ch, 0, 0.0, 0.0, 0.0, qoc);
}


//----------------------------------------
//
//  Add photons to a photon bunch with physical
//  parameters about the photon shape
//
//----------------------------------------
void ph_bunch::add_photons(int N, int ch, int P, double t, double f, double w, qocircuit *qoc){
//  int    N;        // Number of photons to be added
//  int   ch;        // Chanel where the photons are added
//  int    P;        // Polarization of the photons
//  double t;        //Photon emission time
//  double f;        // Photon emission frequency
//  double w;        // Photon emission  width or decay time depending on the packet shape model
//  qocircuit* qoc;  // Circuit to which we will emit the photons
//  Variables
    int    T;        // Packet number
    int    add;      // Add new packer 0=No/1=Yes. This way we check repetitions.
    int    cket;     // Current ket. We may add various photons in different instructions to the same ket
    int   *occ;      // Level occupation
    state *aux;      // Auxiliary state
    hterm  in_term;  // Human readable input term
    state *newbunch; // New bunch of photons with the new photon information updated with respect the old one
// Auxiliary index
    int    i;        // Aux index


    // Check if the packet described for this entry already exists
    add=1;
    i=0;
    while((i<npack)&&(add=1)){
        if((pack_list(1,i)==t)&&(pack_list(2,i)==f)&&(pack_list(3,i)==w)) add=0;
        i++;
    }

    // Update the wave-packet information
    if(add==1){
        // If not exists add it
        T=npack;
        pack_list(0,npack)=(double) npack;
        pack_list(1,npack)=         t;
        pack_list(2,npack)=         f;
        pack_list(3,npack)=         w;
        npack=npack+1;
    }else{
        // If exist just take the index
        T=i-1;
    }


    // Update the input state information
    // Compute the new state
    cket=bunch->nket-1;
    in_term.resize(4,1);
    in_term << ch,
               P,
               T,
               N;
    aux=new state(bunch->nlevel,1);
    aux->add_term(1.0,in_term,qoc);
    occ=new int[bunch->nlevel]();

    for(i=0;i<bunch->nlevel;i++){
        occ[i]=bunch->ket[cket][i]+aux->ket[0][i];
    }

    // Store the new resulting state
    newbunch=new state(bunch->nlevel,bunch->maxket);
    for(i=0;i<cket;i++){
        newbunch->add_term(1.0,bunch->ket[i]);
    }
    newbunch->add_term(1.0,occ);

    //Update internally the bunch
    delete bunch;
    bunch=newbunch;
}


//----------------------------------------
//
// Amplitude of the current alternative/ket
//
//----------------------------------------
void ph_bunch::weight(cmplx A){
//  cmplx A;        // Amplitude


    bunch->ampl[bunch->nket-1]=A;
}


//----------------------------------------
//
//  It creates a new ket that represent a
//  superposition of states
//
//----------------------------------------
void ph_bunch::alternative(){
    int *occ;   // Level occupation


    occ=new int[bunch->nlevel]();
    bunch->add_term(1.0,occ);

    delete occ;
}


//----------------------------------------
//
//  Send the photons to the circuit.
//  Configure th emitter in the circuit.
//  It uses default values for the photons
//  shape model
//
//----------------------------------------
void ph_bunch::send2circuit(qocircuit *qoc){
//  qocircuit *qoc;     // Circuit to which the bunch is referred.


    send2circuit('G',0,qoc);
}


//----------------------------------------
//
//  Send the photons to the circuit specifying
//  the photons shape model. Configure the emitter.
//
//----------------------------------------
void ph_bunch::send2circuit(char ckind, int rand,qocircuit *qoc){
//  char ckind;         // Photon shape model 'G': Gaussian/ 'E': Exponential
//  int  rand;          // Emission time randomization model
//  qocircuit *qoc;     // Circuit to which the bunch is referred.


    qoc->pack_list=pack_list;
    qoc->npack=npack;
    if(qoc->ns>1) qoc->emitter(ckind,rand); // If we only have one possible packer emitter do not perform Gram-Schmidt well
    bunch->normalize();
}


//----------------------------------------
//
//  Returns the internal photon state.
//
//----------------------------------------
state *ph_bunch::bunch_state(){


    return bunch;
}


//----------------------------------------
//
//  Prints internal state
//
//----------------------------------------
void ph_bunch::prnt_state(){


    bunch->prnt_state(nullptr,false,0);
}


//----------------------------------------
//
//  Prints internal state
//
//----------------------------------------

void  ph_bunch::prnt_state(qocircuit *qoc, int column){
//  qocircuit *qoc;     // Circuit to which the ket is referred. (Necessary to print in human readable form)
//  int column          // Do you want the state print in line or with kets in columns? 0=No/1=Yes


    bunch->prnt_state(qoc, false, column);
}


//----------------------------------------
//
//  Prints internal state
//
//----------------------------------------
void  ph_bunch::prnt_state(qocircuit *qoc, bool loss,int column){
//  qocircuit *qoc;     // Circuit to which the ket is referred. (Necessary to print in human readable form)
//  bool loss;          // Print loss channels in different color? False=No/True=Yes
//  int column          // Do you want the state print in line or with kets in columns? 0=No/1=Yes


    bunch->prnt_state(qoc, loss,column);
}
