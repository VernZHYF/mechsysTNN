/************************************************************************
 * MechSys - Open Library for Mechanical Systems                        *
 * Copyright (C) 2009 Sergio Galindo                                    *
 *                                                                      *
 * This program is free software: you can redistribute it and/or modify *
 * it under the terms of the GNU General Public License as published by *
 * the Free Software Foundation, either version 3 of the License, or    *
 * any later version.                                                   *
 *                                                                      *
 * This program is distributed in the hope that it will be useful,      *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of       *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the         *
 * GNU General Public License for more details.                         *
 *                                                                      *
 * You should have received a copy of the GNU General Public License    *
 * along with this program. If not, see <http://www.gnu.org/licenses/>  *
 ************************************************************************/
// Colloid transport

//STD
#include<iostream>

// MechSys
#include <mechsys/lbm/Domain.h>
#include <mechsys/util/maps.h>
#include <mechsys/util/util.h>

struct UserData
{
    Array<Cell *> Left; 
    Array<Cell *> Right;
    size_t         alt;
    size_t          Np;
    double        vmax; 
    double          Kn; 
    double        xlim;
    double          rc;
    double          vp;
    Vec3_t        Xmin; 
    Vec3_t        Xmax; 
    Vec3_t           g;
    std::ofstream      oss_ss;       ///< file for stress strain data
};

void Setup (LBM::Domain & dom, void * UD)
{
    UserData & dat = (*static_cast<UserData *>(UD));
    //for (size_t i=0;i<dat.Left.Size();i++)
    //{
        //dat.Left [i]->Initialize(dat.Left [i]->RhoBC,dat.Left [i]->VelBC);
        //dat.Right[i]->Initialize(dat.Right[i]->RhoBC,dat.Right[i]->VelBC);
    //}

	// Cells with prescribed velocity
    #pragma omp parallel for schedule(static) num_threads(dom.Nproc)
	for (size_t i=0; i<dat.Left.Size(); ++i)
	{
		Cell * c = dat.Left[i];
		if (c->IsSolid) continue;
		double rho = (c->F[0]+c->F[2]+c->F[4] + 2.0*(c->F[3]+c->F[6]+c->F[7]))/(1.0-c->VelBC(0));
		c->F[1] = c->F[3] + (2.0/3.0)*rho*c->VelBC(0);
		c->F[5] = c->F[7] + (1.0/6.0)*rho*c->VelBC(0) + 0.5*rho*c->VelBC(1) - 0.5*(c->F[2]-c->F[4]);
		c->F[8] = c->F[6] + (1.0/6.0)*rho*c->VelBC(0) - 0.5*rho*c->VelBC(1) + 0.5*(c->F[2]-c->F[4]);
        c->Rho = c->VelDen(c->Vel);
	}

	// Cells with prescribed density 
    #pragma omp parallel for schedule(static) num_threads(dom.Nproc)
	for (size_t i=0; i<dat.Right.Size(); ++i)
	{
		Cell * c = dat.Right[i];
		if (c->IsSolid) continue;
		double vx = -1.0 + (c->F[0]+c->F[2]+c->F[4] + 2.0*(c->F[1]+c->F[5]+c->F[8]))/c->RhoBC;
		c->F[3] = c->F[1] - (2.0/3.0)*c->RhoBC*vx; 
		c->F[7] = c->F[5] - (1.0/6.0)*c->RhoBC*vx + 0.5*(c->F[2]-c->F[4]);
		c->F[6] = c->F[8] - (1.0/6.0)*c->RhoBC*vx - 0.5*(c->F[2]-c->F[4]);
        c->Rho = c->VelDen(c->Vel);
	}

    dat.Np = 0;
    #pragma omp parallel for schedule(static) num_threads(dom.Nproc)
    for (size_t i=0;i<dom.Disks.Size();i++)
    {
        if (dom.Disks[i]->X(0) > dat.Xmax(0)) dat.Np++;
        dom.Disks[i]->Ff = 0.0,0.0,0.0;
        double delta;
        delta =   dat.Xmin(0) - dom.Disks[i]->X(0) + dom.Disks[i]->R;
        if (delta > 0.0)  dom.Disks[i]->Ff(0) += dat.Kn*delta;
        delta =   dat.Xmin(1) - dom.Disks[i]->X(1) + dom.Disks[i]->R;
        if (delta > 0.0)  dom.Disks[i]->Ff(1) += dat.Kn*delta;
        delta = - dat.Xmax(1) + dom.Disks[i]->X(1) + dom.Disks[i]->R;
        if (delta > 0.0)  dom.Disks[i]->Ff(1) -= dat.Kn*delta;
    }

    Vec3_t Xmin,Xmax;
    dom.BoundingBox(Xmin,Xmax);
    if (Xmin(0) > dat.xlim + 2.0*dat.rc)
    {
        for (double y=0.05*dat.Xmax(1) + dat.alt*2.0*dat.rc;y<0.95*dat.Xmax(1);y+=4.0*dat.rc)
        {
            dom.AddDisk(0,Vec3_t(dat.xlim, y,0.0),Vec3_t(dat.vp,0.0,0.0),OrthoSys::O,3.0,dat.rc,1.0);
            dom.Disks[dom.Disks.Size()-1]->Kn = dat.Kn;
        }
        dom.ResetContacts();
        dom.ResetDisplacements();
        dat.alt ++;
        dat.alt = (dat.alt%2);
    }
}

void Report(LBM::Domain & dom, void * UD)
{
    UserData & dat = (*static_cast<UserData *>(UD));
    if (dom.idx_out==0)
    {
        String fs;
        fs.Printf("%s_n_particles.res","test04");
        dat.oss_ss.open(fs.CStr(),std::ios::out);
        dat.oss_ss << "Time" << Util::_8s << "Np" <<"\n";
    }

    dat.oss_ss << dom.Time << Util::_8s << dat.Np << "\n";
}

int main(int argc, char **argv) try
{
    size_t Nproc  = 1; 
    size_t nx     = 1000;
    size_t ny     = 1000;
    double nu     = 0.04;
    double dx     = 1.0;
    double dt     = 1.0;    
    double vb     = 0.01;
    double rc     = 5.0;
    double xlim   = 1.1*rc;
    double Kn     = 1.0e2;
    double Tf     = 1.0e5;
    if (argc>=2) Nproc = atoi(argv[1]);
    if (argc>=3) Tf    = atof(argv[2]);

    LBM::Domain Dom(D2Q9, nu, iVec3_t(nx,ny,1), dx, dt);
    UserData dat;
    Dom.UserData = &dat;
    Dom.Alpha    = 2.0*rc;
    Dom.Sc       = 0.0;
    dat.xlim     = xlim;
    dat.Xmin     = 0.0,0.0,0.0;
    dat.Xmax     = nx*dx,ny*dx,0.0;
    dat.rc       = rc;
    dat.vp       = 0.0*vb;
    dat.g        = 0.001,0.0,0.0;
    dat.alt      = 1;
    dat.Np       = 0;
    dat.Kn       = Kn;

    //std::cout << "1" << std::endl;
    //Assigning the left and right cells
    for (size_t i=0;i<ny;i++)
    {
        dat.Left .Push(Dom.Lat[0].GetCell(iVec3_t(0   ,i,0)));
        dat.Right.Push(Dom.Lat[0].GetCell(iVec3_t(nx-1,i,0)));
        double vx = vb; // horizontal velocity
        double vy = 0.0;                          // vertical velocity
		Vec3_t v(vx, vy, 0.0);                    // velocity vector
        dat.Left [i]->VelBC = v;
        dat.Left [i]->RhoBC = 1.0;
        dat.Right[i]->VelBC = v;
        dat.Right[i]->RhoBC = 1.0;

    }

    //Assigning solid boundaries at top and bottom
    for (size_t i=0;i<nx;i++)
    {
        Dom.Lat[0].GetCell(iVec3_t(i,0   ,0))->IsSolid = true;
        Dom.Lat[0].GetCell(iVec3_t(i,ny-1,0))->IsSolid = true;
    }

    //std::cout << "2" << std::endl;
	// Set grains
	Table grains;
	grains.Read("circles.out");
	for (size_t i=0; i<grains["Xc"].Size(); ++i)
	{
		double xc = grains["Xc"][i]*nx+0.05*nx;
		double yc = grains["Yc"][i]*ny;
		double r  = grains["R" ][i]*nx*0.9;
        Dom.AddDisk(0,Vec3_t(xc,yc,0.0),OrthoSys::O,OrthoSys::O,3.0,r,1.0);
        Dom.Disks[Dom.Disks.Size()-1]->FixVeloc();
        Dom.Disks[Dom.Disks.Size()-1]->Kn = Kn;
	}

    //std::cout << "3" << std::endl;
    for (double y=0.05*dat.Xmax(1);y<0.95*dat.Xmax(1);y+=4.0*rc)
    {
        Dom.AddDisk(0,Vec3_t(xlim, y,0.0),Vec3_t(0.1*vb,0.0,0.0),OrthoSys::O,3.0,rc,1.0);
        Dom.Disks[Dom.Disks.Size()-1]->Kn = Kn;
    }
    //Dom.AddDisk(0,Vec3_t(0.3*nx,0.5*ny,0.0),Vec3_t(0.0,0.0,0.0),OrthoSys::O,3.0,10.0*rc,1.0);

    

    double rho0 = 1.0;
    Vec3_t v0(0.0,0.0,0.0);

    //Initializing values
    for (size_t i=0;i<Dom.Lat[0].Ncells;i++)
    {
        Dom.Lat[0].Cells[i]->Initialize(rho0, v0);
    }


    //Solving
    Dom.Solve(Tf,0.01*Tf,Setup,Report,"test04",true,Nproc);
    dat.oss_ss.close();
 
}
MECHSYS_CATCH

