/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2012 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "StreamParcel.H"
#include "IOstreams.H"
#include "IOField.H"
#include "Cloud.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

template<class ParcelType>
Foam::string Foam::StreamParcel<ParcelType>::propertyList_ =
    Foam::StreamParcel<ParcelType>::propertyList();

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ParcelType>
Foam::StreamParcel<ParcelType>::StreamParcel
(
    const polyMesh& mesh,
    Istream& is,
    bool readFields
)
:
    ParcelType(mesh, is, readFields),
    active_(false),
    typeId_(0),
    nParticle_(0.0),
    d_(0.0),
    dTarget_(0.0),
    U_(vector::zero),
    f_(vector::zero),
    angularMomentum_(vector::zero),
    torque_(vector::zero),
    rho_(0.0),
    age_(0.0),
    cV_(0.0),
    //- modified by Patrick
    //Mfp_(0.0),
    //Tp_(300),
    //YTTIPp_(0.0), 
    //RRp_(0.0), 
    N_(0.0),
    A_(0.0),
    V_(0.0),
    a_(0.0),
    v_(0.0),    
    dp_(0.0),
    rc_(0.0),
    cn_(0.0),
    Dp_(0.0),
    Kn_(0.0),
    Cun_(0.0),	
    l_(0.0),
    gp_(0.0),
    beta_(0.0),
    tau_(GREAT),       
    tTurb_(0.0),
    UTurb_(vector::zero),
    rhoc_(0.0),
    Uc_(vector::zero),
    muc_(0.0),
    //Mfc_(0.0),
    //epsilonc_(0.0),
    //kc_(0.0),
    //YTTIPc_(0.0)
    //- added by Patrick
    Tc_(300.0),
    concSourcec_(0.0),
    gradTc_(vector::zero)	
{
    if (readFields)
    {
        if (is.format() == IOstream::ASCII)
        {
            active_ = readBool(is);
            typeId_ = readLabel(is);
            nParticle_ = readScalar(is);
            d_ = readScalar(is);
            dTarget_ = readScalar(is);
            is >> U_;
            is >> f_;
            is >> angularMomentum_;
            is >> torque_;
            rho_ = readScalar(is);
            age_ = readScalar(is);
	    cV_ = readScalar(is);
	    //- modified by Patrick
            //Mfp_ = readScalar(is); 
            //Tp_ = readScalar(is);
            //YTTIPp_ = readScalar(is);	
            //RRp_ = readScalar(is);
            N_= readScalar(is);
	    A_= readScalar(is);
	    V_= readScalar(is);
            a_= readScalar(is);
            v_= readScalar(is);    
            dp_= readScalar(is);
            rc_= readScalar(is);
            cn_= readScalar(is);
            Dp_= readScalar(is);
            Kn_= readScalar(is);
            Cun_= readScalar(is);	    
            l_= readScalar(is);
            gp_= readScalar(is);
            beta_= readScalar(is);
            tau_= readScalar(is);     	    
            tTurb_ = readScalar(is);
            is >> UTurb_;
	    Tc_ = readScalar(is);
	    concSourcec_ = readScalar(is);
	    is >> gradTc_;
        }
        else
        {
            is.read
            (
                reinterpret_cast<char*>(&active_),
                sizeof(active_)
              + sizeof(typeId_)
              + sizeof(nParticle_)
              + sizeof(d_)
              + sizeof(dTarget_)
              + sizeof(U_)
              + sizeof(f_)
              + sizeof(angularMomentum_)
              + sizeof(torque_)
              + sizeof(rho_)
              + sizeof(age_)
	      + sizeof(cV_)
	      //- modified by Patrick
              //+ sizeof(Mfp_)  
              //+ sizeof(Tp_)
              //+ sizeof(YTTIPp_)
              //+ sizeof(RRp_)
              + sizeof(N_)
              + sizeof(A_)
              + sizeof(V_)
              + sizeof(a_)
              + sizeof(v_)    
              + sizeof(dp_)
              + sizeof(rc_)
              + sizeof(cn_)
              + sizeof(Dp_)
	      + sizeof(Kn_)
	      + sizeof(Cun_)
              + sizeof(l_)
              + sizeof(gp_)
              + sizeof(beta_)
              + sizeof(tau_)  	      
              + sizeof(tTurb_)
              + sizeof(UTurb_)
	      + sizeof(Tc_)
	      + sizeof(concSourcec_)
	      + sizeof(gradTc_)
            );
        }
    }

    // Check state of Istream
    is.check
    (
        "StreamParcel<ParcelType>::StreamParcel"
        "(const polyMesh&, Istream&, bool)"
    );
}


template<class ParcelType>
template<class CloudType>
void Foam::StreamParcel<ParcelType>::readFields(CloudType& c)
{
    if (!c.size())
    {
        return;
    }

    ParcelType::readFields(c);

    IOField<label> active(c.fieldIOobject("active", IOobject::MUST_READ));
    c.checkFieldIOobject(c, active);

    IOField<label> typeId(c.fieldIOobject("typeId", IOobject::MUST_READ));
    c.checkFieldIOobject(c, typeId);

    IOField<scalar>
        nParticle(c.fieldIOobject("nParticle", IOobject::MUST_READ));
    c.checkFieldIOobject(c, nParticle);

    IOField<scalar> d(c.fieldIOobject("d", IOobject::MUST_READ));
    c.checkFieldIOobject(c, d);

    IOField<scalar> dTarget(c.fieldIOobject("dTarget", IOobject::MUST_READ));
    c.checkFieldIOobject(c, dTarget);

    IOField<vector> U(c.fieldIOobject("U", IOobject::MUST_READ));
    c.checkFieldIOobject(c, U);

    IOField<vector> f(c.fieldIOobject("f", IOobject::MUST_READ));
    c.checkFieldIOobject(c, f);

    IOField<vector> angularMomentum
    (
        c.fieldIOobject("angularMomentum", IOobject::MUST_READ)
    );
    c.checkFieldIOobject(c, angularMomentum);

    IOField<vector> torque(c.fieldIOobject("torque", IOobject::MUST_READ));
    c.checkFieldIOobject(c, torque);

    IOField<scalar> rho(c.fieldIOobject("rho", IOobject::MUST_READ));
    c.checkFieldIOobject(c, rho);

    IOField<scalar> age(c.fieldIOobject("age", IOobject::MUST_READ));
    c.checkFieldIOobject(c, age);
    
 
    IOField<scalar> cV(c.fieldIOobject("cV", IOobject::MUST_READ));
    c.checkFieldIOobject(c, cV);   

    //- modified by Patrick

    //IOField<scalar> Mfp(c.fieldIOobject("Mfp", IOobject::MUST_READ));
    //c.checkFieldIOobject(c, Mfp);   

    //IOField<scalar> Tp(c.fieldIOobject("Tp", IOobject::MUST_READ));
    //c.checkFieldIOobject(c, Tp);     

    //IOField<scalar> YTTIPp(c.fieldIOobject("YTTIPp", IOobject::MUST_READ));
    //c.checkFieldIOobject(c, YTTIPp);  

    //IOField<scalar> RRp(c.fieldIOobject("RRp", IOobject::MUST_READ));
    //c.checkFieldIOobject(c, RRp);  

    IOField<scalar> N(c.fieldIOobject("N", IOobject::MUST_READ));
    c.checkFieldIOobject(c, N);  
    

    IOField<scalar> A(c.fieldIOobject("A", IOobject::MUST_READ));
    c.checkFieldIOobject(c, A);  

    IOField<scalar> V(c.fieldIOobject("V", IOobject::MUST_READ));
    c.checkFieldIOobject(c, V);      

    IOField<scalar> a(c.fieldIOobject("a", IOobject::MUST_READ));
    c.checkFieldIOobject(c, a);      

    IOField<scalar> v(c.fieldIOobject("v", IOobject::MUST_READ));
    c.checkFieldIOobject(c, v);      

    IOField<scalar> dp(c.fieldIOobject("dp", IOobject::MUST_READ));
    c.checkFieldIOobject(c, dp);      

    IOField<scalar> rc(c.fieldIOobject("rc", IOobject::MUST_READ));
    c.checkFieldIOobject(c, rc);      

    IOField<scalar> cn(c.fieldIOobject("cn", IOobject::MUST_READ));
    c.checkFieldIOobject(c, cn); 

    IOField<scalar> Dp(c.fieldIOobject("Dp", IOobject::MUST_READ));
    c.checkFieldIOobject(c, Dp); 
   
    IOField<scalar> Kn(c.fieldIOobject("Kn", IOobject::MUST_READ));
    c.checkFieldIOobject(c, Kn); 

    IOField<scalar> Cun(c.fieldIOobject("Cun", IOobject::MUST_READ));
    c.checkFieldIOobject(c, Cun); 

    IOField<scalar> l(c.fieldIOobject("l", IOobject::MUST_READ));
    c.checkFieldIOobject(c, l);      

    IOField<scalar> gp(c.fieldIOobject("gp", IOobject::MUST_READ));
    c.checkFieldIOobject(c, gp);    

    IOField<scalar> beta(c.fieldIOobject("beta", IOobject::MUST_READ));
    c.checkFieldIOobject(c, beta);       
 
    IOField<scalar> tau(c.fieldIOobject("tau", IOobject::MUST_READ));
    c.checkFieldIOobject(c, tau);      
        
    IOField<scalar> tTurb(c.fieldIOobject("tTurb", IOobject::MUST_READ));
    c.checkFieldIOobject(c, tTurb);

    IOField<vector> UTurb(c.fieldIOobject("UTurb", IOobject::MUST_READ));
    c.checkFieldIOobject(c, UTurb);

    IOField<scalar> Tc(c.fieldIOobject("Tc", IOobject::MUST_READ));
    c.checkFieldIOobject(c, Tc);    

    IOField<scalar> concSourcec(c.fieldIOobject("concSourcec", IOobject::MUST_READ));
    c.checkFieldIOobject(c, concSourcec);      

    IOField<vector> gradTc(c.fieldIOobject("gradTc", IOobject::MUST_READ));
    c.checkFieldIOobject(c, gradTc);

    label i = 0;

    forAllIter(typename CloudType, c, iter)
    {
        StreamParcel<ParcelType>& p = iter();

        p.active_ = active[i];
        p.typeId_ = typeId[i];
        p.nParticle_ = nParticle[i];
        p.d_ = d[i];
        p.dTarget_ = dTarget[i];
        p.U_ = U[i];
        p.f_ = f[i];
        p.angularMomentum_ = angularMomentum[i];
        p.rho_ = rho[i];
        p.age_ = age[i];
	p.cV_ = cV[i];
        //p.Mfp_ = Mfp[i];
        //p.Tp_ = Tp[i];
        //p.YTTIPp_ = YTTIPp[i];
        //p.RRp_ = RRp[i];
	p.N_ = N[i]; 
	p.A_ = A[i]; 
	p.V_ = V[i]; 
	p.a_ = a[i]; 
	p.v_ = v[i]; 
	p.dp_ = dp[i]; 
	p.rc_ = rc[i]; 
	p.cn_ = cn[i]; 
	p.Dp_ = Dp[i];
        p.Kn_ = Kn[i];
        p.Cun_ = Cun[i];	
	p.l_ = l[i]; 
	p.gp_ = gp[i]; 
	p.beta_ = beta[i]; 
	p.tau_ = tau[i]; 
        p.tTurb_ = tTurb[i];
        p.UTurb_ = UTurb[i];
	p.Tc_ = Tc[i];
	p.concSourcec_ = concSourcec[i];
	p.gradTc_ = gradTc[i];

        i++;
    }
}


template<class ParcelType>
template<class CloudType>
void Foam::StreamParcel<ParcelType>::writeFields(const CloudType& c)
{
    ParcelType::writeFields(c);

    label np =  c.size();

    IOField<label> active(c.fieldIOobject("active", IOobject::NO_READ), np);
    IOField<label> typeId(c.fieldIOobject("typeId", IOobject::NO_READ), np);
    IOField<scalar> nParticle
    (
        c.fieldIOobject("nParticle", IOobject::NO_READ),
        np
    );
    IOField<scalar> d(c.fieldIOobject("d", IOobject::NO_READ), np);
    IOField<scalar> dTarget(c.fieldIOobject("dTarget", IOobject::NO_READ), np);
    IOField<vector> U(c.fieldIOobject("U", IOobject::NO_READ), np);
    IOField<vector> f(c.fieldIOobject("f", IOobject::NO_READ), np);
    IOField<vector> angularMomentum
    (
        c.fieldIOobject("angularMomentum", IOobject::NO_READ),
        np
    );
    IOField<vector> torque(c.fieldIOobject("torque", IOobject::NO_READ), np);
    IOField<scalar> rho(c.fieldIOobject("rho", IOobject::NO_READ), np);
    IOField<scalar> age(c.fieldIOobject("age", IOobject::NO_READ), np);
    IOField<scalar> cV(c.fieldIOobject("cV", IOobject::NO_READ), np);  
    //- modified by Patrick
    //IOField<scalar> Mfp(c.fieldIOobject("Mfp", IOobject::NO_READ), np); 
    //IOField<scalar> Tp(c.fieldIOobject("Tp", IOobject::NO_READ), np);  
    //IOField<scalar> YTTIPp(c.fieldIOobject("YTTIPp", IOobject::NO_READ), np);
    //IOField<scalar> RRp(c.fieldIOobject("RRp", IOobject::NO_READ), np);   
    IOField<scalar> N(c.fieldIOobject("N", IOobject::NO_READ), np);  
    IOField<scalar> A(c.fieldIOobject("A", IOobject::NO_READ), np);  
    IOField<scalar> V(c.fieldIOobject("V", IOobject::NO_READ), np);  
    IOField<scalar> a(c.fieldIOobject("a", IOobject::NO_READ), np); 
    IOField<scalar> v(c.fieldIOobject("v", IOobject::NO_READ), np);      
    IOField<scalar> dp(c.fieldIOobject("dp", IOobject::NO_READ), np); 
    IOField<scalar> rc(c.fieldIOobject("rc", IOobject::NO_READ), np);
    IOField<scalar> cn(c.fieldIOobject("cn", IOobject::NO_READ), np); 
    IOField<scalar> Dp(c.fieldIOobject("Dp", IOobject::NO_READ), np);
    IOField<scalar> Kn(c.fieldIOobject("Kn", IOobject::NO_READ), np); 
    IOField<scalar> Cun(c.fieldIOobject("Cun", IOobject::NO_READ), np);     
    IOField<scalar> l(c.fieldIOobject("l", IOobject::NO_READ), np);  
    IOField<scalar> gp(c.fieldIOobject("gp", IOobject::NO_READ), np);  
    IOField<scalar> beta(c.fieldIOobject("beta", IOobject::NO_READ), np);  
    IOField<scalar> tau(c.fieldIOobject("tau", IOobject::NO_READ), np);  
    IOField<scalar> tTurb(c.fieldIOobject("tTurb", IOobject::NO_READ), np);
    IOField<vector> UTurb(c.fieldIOobject("UTurb", IOobject::NO_READ), np);
    IOField<scalar> Tc(c.fieldIOobject("Tc", IOobject::NO_READ), np);  
    IOField<scalar> concSourcec(c.fieldIOobject("concSourcec", IOobject::NO_READ), np);
    IOField<vector> gradTc(c.fieldIOobject("gradTc", IOobject::NO_READ), np);

    label i = 0;

    forAllConstIter(typename CloudType, c, iter)
    {
        const StreamParcel<ParcelType>& p = iter();

        active[i] = p.active();
        typeId[i] = p.typeId();
        nParticle[i] = p.nParticle();
        d[i] = p.d();
        dTarget[i] = p.dTarget();
        U[i] = p.U();
        f[i] = p.f();
        angularMomentum[i] = p.angularMomentum();
        torque[i] = p.torque();
        rho[i] = p.rho();
        age[i] = p.age();
	cV[i] = p.cV();
	//- modified by Patrick
        //Mfp[i] = p.Mfp();
        //Tp[i] = p.Tp();	
        //YTTIPp[i] = p.YTTIPp();	
        //RRp[i] = p.RRp();
	N[i]=p.N(); 
	A[i]=p.A(); 
	V[i]=p.V(); 
	a[i]=p.a(); 
	v[i]=p.v(); 
	dp[i]=p.dp(); 
	rc[i]=p.rc(); 
	cn[i]=p.cn();
	Dp[i]=p.Dp();
	Kn[i]=p.Kn();
	Cun[i]=p.Cun();
	l[i]=p.l(); 
	gp[i]=p.gp(); 
	 beta[i]=p.beta(); 
	tau[i]=p.tau(); 	
        tTurb[i] = p.tTurb();
        UTurb[i] = p.UTurb();
	Tc[i] = p.Tc();
	concSourcec[i] = p.concSourcec();
	gradTc[i] = p.gradTc();

        i++;
    }

    active.write();
    typeId.write();
    nParticle.write();
    d.write();
    dTarget.write();
    U.write();
    f.write();
    angularMomentum.write();
    torque.write();
    rho.write();
    age.write();
    cV.write();
    //- modified by Patrick
    //Mfp.write();
    //Tp.write();
    //YTTIPp.write();
    //RRp.write();    
    N.write();
    A.write(); 
    V.write(); 
    a.write();
    v.write();
    dp.write(); 
    rc.write(); 
    cn.write();
    Dp.write();
    Kn.write();
    Cun.write();
    l.write();
    gp.write(); 
    beta.write();
    tau.write(); 	        
    tTurb.write();
    UTurb.write();
    Tc.write();
    concSourcec.write();
    gradTc.write();
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

template<class ParcelType>
Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const StreamParcel<ParcelType>& p
)
{
    if (os.format() == IOstream::ASCII)
    {
        os  << static_cast<const ParcelType&>(p)
            << token::SPACE << p.active()
            << token::SPACE << p.typeId()
            << token::SPACE << p.nParticle()
            << token::SPACE << p.d()
            << token::SPACE << p.dTarget()
            << token::SPACE << p.U()
            << token::SPACE << p.f()
            << token::SPACE << p.angularMomentum()
            << token::SPACE << p.torque()
            << token::SPACE << p.rho()
            << token::SPACE << p.age()
	    << token::SPACE << p.cV()
	    //- modified by Patrick
            //<< token::SPACE << p.Mfp()
            //<< token::SPACE << p.Tp()
            //<< token::SPACE << p.YTTIPp()
            //<< token::SPACE << p.RRp()
            << token::SPACE << p.N()
            << token::SPACE << p.A()	    
            << token::SPACE << p.V()	    
            << token::SPACE << p.a()
            << token::SPACE << p.v()	    
            << token::SPACE << p.dp()	    
            << token::SPACE << p.rc()	
            << token::SPACE << p.cn()
            << token::SPACE << p.Dp()
	    << token::SPACE << p.Kn()
	    << token::SPACE << p.Cun()
            << token::SPACE << p.l()
            << token::SPACE << p.gp()	    
            << token::SPACE << p.beta()	
            << token::SPACE << p.tau()		    
            << token::SPACE << p.tTurb()	    
            << token::SPACE << p.UTurb()
	    << token::SPACE << p.Tc()
	    << token::SPACE << p.concSourcec()
	    << token::SPACE << p.gradTc();
    }
    else
    {
        os  << static_cast<const ParcelType&>(p);
        os.write
        (
            reinterpret_cast<const char*>(&p.active_),
            sizeof(p.active())
          + sizeof(p.typeId())
          + sizeof(p.nParticle())
          + sizeof(p.d())
          + sizeof(p.dTarget())
          + sizeof(p.U())
          + sizeof(p.f())
          + sizeof(p.angularMomentum())
          + sizeof(p.torque())
          + sizeof(p.rho())
          + sizeof(p.age())
	  + sizeof(p.cV())
	  //- modified by Patrick
          //+ sizeof(p.Mfp()) 
          //+ sizeof(p.Tp())
          //+ sizeof(p.YTTIPp())
          //+ sizeof(p.RRp()) 
          + sizeof(p.N())
          + sizeof(p.A())
          + sizeof(p.V())
          + sizeof(p.a())
          + sizeof(p.v())    
          + sizeof(p.dp())
          + sizeof(p.rc())
          + sizeof(p.cn())
          + sizeof(p.Dp())
	  + sizeof(p.Kn())
	  + sizeof(p.Cun())
          + sizeof(p.l())
          + sizeof(p.gp())
          + sizeof(p.beta())
          + sizeof(p.tau()) 
	  + sizeof(p.tTurb())
          + sizeof(p.UTurb())
	  + sizeof(p.Tc())
	  + sizeof(p.concSourcec())
	  + sizeof(p.gradTc())
        );
    }

    // Check state of Ostream
    os.check
    (
        "Ostream& operator<<(Ostream&, const StreamParcel<ParcelType>&)"
    );

    return os;
}


// ************************************************************************* //
