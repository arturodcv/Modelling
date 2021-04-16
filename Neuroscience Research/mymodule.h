/*
 *  mymodule.h
 *
 *  This file is part of NEST.
 *
 *  Copyright (C) 2004 The NEST Initiative
 *
 *  NEST is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  NEST is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with NEST.  If not, see <http://www.gnu.org/licenses/>.
 *
 */
 
 
 //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
 //
 // Module created by Arturo del Cerro in order to replicate the work shown in 
 //'A Neurodynamical Model of Brightness Induction in V1' by Olivier Penacchio, 
 // Xavier Otazu and Laura Dempere-Marco.
 //
 // This module is GNU General Public License.
 //
 //////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#ifndef MYMODULE_H
#define MYMODULE_H

// Includes from sli:
#include "topologymodule.h"
#include "topology_parameter.h"
#include "slifunction.h"
#include "slimodule.h"
#include <math.h>

// Put your stuff into your own namespace.
namespace mynest
{

/**
 * Class defining your model.
 * @note For each model, you must define one such class, with a unique name.
 */
class MyModule : public SLIModule
{
public:
  // Interface functions ------------------------------------------

  /**
   * @note The constructor registers the module with the dynamic loader.
   *       Initialization proper is performed by the init() method.
   */
  MyModule();

  /**
   * @note The destructor does not do much in modules.
   */
  ~MyModule();

  /**
   * Initialize module.
   * @param SLIInterpreter* SLI interpreter
   */
  void init( SLIInterpreter* );

  /**
   * Return the name of your model.
   */
  const std::string name( void ) const;

  /**
   * Return the name of a sli file to execute when mymodule is loaded.
   * This mechanism can be used to define SLI commands associated with your
   * module, in particular, set up type tries for functions you have defined.
   */
  const std::string commandstring( void ) const;
};

class RadialParameter : public nest::TopologyParameter
{
public:
  RadialParameter()
    : TopologyParameter()
  {
  }

  RadialParameter( double cutoff )
    : TopologyParameter( cutoff )
  {
  }

  RadialParameter( const DictionaryDatum& d )
    : TopologyParameter( d )
  {
  }

  virtual double raw_value( double ) const = 0;

  double
  raw_value( const nest::Position< 2 >& p, librandom::RngPtr& ) const
  {
    return raw_value( p.length() );
  }
  double
  raw_value( const nest::Position< 3 >& p, librandom::RngPtr& ) const
  {
    return raw_value( p.length() );
  }
};

class PlosOne_J: public nest::TopologyParameter
{
public:
  PlosOne_J(const DictionaryDatum& d):
    TopologyParameter(d),
    rows_(1.0),
    kappa(1.0),
    orientation_i(0.0),
    orientation_j(0.0)
    {
      updateValue<double>(d, "rows", rows_);
      updateValue<double>(d, "kappa", kappa);
      updateValue<double>(d, "orientation_i", orientation_i);
      updateValue<double>(d, "orientation_j", orientation_j);	
    }

  double raw_value(const nest::Position<2>& pos, librandom::RngPtr&) const
    { 
	double theta_1,theta_2,theta_aux,beta,dist,exp_,angle,orientation_1,orientation_2;
        double pi = M_PI; 
	double euler = exp(1.0);

	angle = atan2( pos[1], pos[0]);



        orientation_1 = orientation_i;
	orientation_2 = orientation_j;


	if (orientation_1 < 0.0){orientation_1 = orientation_1 + pi;}
	if (orientation_2 < 0.0){orientation_2 = orientation_2 + pi;}

 
	theta_1 = fabs(orientation_1 - angle);
  theta_2 = fabs(orientation_2 - angle);

	if (fabs(theta_1) > pi/2){ theta_1 = -(theta_1 - pi) ;}
        if (fabs(theta_2) > pi/2){ theta_2 = -(theta_2 - pi) ;}  

	if (theta_1 <= - pi/2) { theta_1 = -theta_1  ;}
        if (theta_2 <= - pi/2) { theta_2 = -theta_2  ;}

        
    	if (fabs(theta_1) > fabs(theta_2)){
        	theta_aux = theta_2;
        	theta_2 = theta_1;
        	theta_1 = theta_aux;
	}

	
	beta = 2 * ( fabs(theta_1) + sin ( fabs ( theta_1 + theta_2 )));


	dist = sqrt( pow(pos[1],2) + pow(pos[0],2) ) * rows_;

	

	if ( (dist > 0.0 && dist <= 10.0 && beta < pi/2.69) || (
	   ( (dist > 0.0 && dist <= 10.0 && beta < pi/1.1) && (fabs(theta_1) < pi/5.9)) && (fabs(theta_2) < pi/5.9) )){

         	exp_ = - pow( beta / dist , 2)  - 2 * pow ( beta / dist,7) - pow( dist , 2) / 90;
		
		return  kappa * pow( euler , exp_);}

	else { return 0.0;}
    }

  nest::TopologyParameter * clone() const
    { return new PlosOne_J(*this); }

private:
  double rows_,kappa, orientation_i, orientation_j;
};




class PlosOne_W: public nest::TopologyParameter
{
public:
  PlosOne_W(const DictionaryDatum& d):
    TopologyParameter(d),
    rows_(1.0),
    orientation_i(0.0),
    orientation_j(0.0),
    kappa(1.0)
    {
      updateValue<double>(d, "rows", rows_);
      updateValue<double>(d, "orientation_i", orientation_i);
      updateValue<double>(d, "orientation_j", orientation_j);
      updateValue<double>(d, "kappa", kappa);
    }


  double raw_value(const nest::Position<2>& pos,
                   librandom::RngPtr&) const
    { 
	double theta_1,theta_2,theta_aux,beta,dist,exp_1,exp_2,angle,orientation_1,orientation_2;
        double pi = M_PI;
	double euler = exp(1.0);


	angle = atan2( pos[1], pos[0]);


        orientation_1 = orientation_i;
	orientation_2 = orientation_j;
	
	if (orientation_1 < 0.0){orientation_1 = orientation_1 + pi;}
	if (orientation_2 < 0.0){orientation_2 = orientation_2 + pi;}

	double theta_diff = fabs(orientation_1 - orientation_2);
	if (orientation_1 > pi/2 ){
		if (fabs(orientation_1 - orientation_2 - pi) < theta_diff){ theta_diff = fabs(orientation_1 - orientation_2 - pi);}
	}
	if (orientation_2 > pi/2 ){
		if (fabs(orientation_1 - orientation_2 + pi) < theta_diff){ theta_diff = fabs(orientation_1 - orientation_2 + pi);}
	}

	if (theta_diff >= pi){theta_diff = theta_diff - pi;}
	
	theta_1 = fabs(orientation_1 - angle);
        theta_2 = fabs(orientation_2 - angle);


	if (fabs(theta_1) > pi/2){ theta_1 = -(theta_1 - pi) ;}
        if (fabs(theta_2) > pi/2){ theta_2 = -(theta_2 - pi) ;}
	if (theta_1 <= - pi/2) { theta_1 = -theta_1  ;}
        if (theta_2 <= - pi/2) { theta_2 = -theta_2  ;}
        
    	if (fabs(theta_1) > fabs(theta_2)){
        	theta_aux = theta_2;
        	theta_2 = theta_1;
        	theta_1 = theta_aux;
	}
	
	beta = 2*(fabs(theta_1) + sin(fabs(theta_1+theta_2)));

	dist = sqrt(pow(pos[1],2) + pow(pos[0],2)) * rows_;

	if ((dist < 0.1 || dist >= 10.0 || beta < pi / 1.1) || fabs(theta_diff) >= pi/3 || fabs(theta_1) < pi / 11.999){
		return 0.0;}
	else {
		exp_1 = -0.4 * pow( beta / dist , 1.5);
        	exp_2 = -pow( theta_diff / (pi / 4) , 1.5);
	
        	return kappa * ( 1 - pow( euler , exp_1) ) * pow( euler , exp_2);
	}
    }

  nest::TopologyParameter * clone() const
    { return new PlosOne_W(*this); }

private:
  double rows_, orientation_i, orientation_j, kappa;
};



} // namespace mynest

#endif




