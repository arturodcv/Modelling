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
 // Module created by Arturo del Cerro.
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

  const double euler = exp(1.0);
  const double PI_11_999 = M_PI * 0.08333; // 1/11.999 = 0.08333
  const double PI_1_35 = M_PI * 0.7407;
  const double PI_3 = M_PI * 0.3333;
  const double M_4_PI = M_2_PI * 2;
  const double PI_2_69 = M_PI * 0.3717; 
  const double PI_5_9  = M_PI * 0.1694;

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
    kappa(1.0),
    orientation_i(0.0),
    orientation_j(0.0),
    rescale(1.0)
    {
      updateValue<double>(d, "kappa", kappa);
      updateValue<double>(d, "orientation_i", orientation_i);
      updateValue<double>(d, "orientation_j", orientation_j);	
      updateValue<double>(d, "rescale", rescale);
    }

  double raw_value(const nest::Position<2>& pos, librandom::RngPtr&) const
    { 
	double theta_1,theta_2,theta_aux,beta,dist,exp_,angle;

	angle = atan( pos[1] / pos[0]);

	theta_1 = angle - orientation_i; theta_2 = angle - orientation_j;

  if (theta_1 < - M_PI_2){theta_1 = theta_1 + M_PI;}
  if (theta_1 > M_PI_2){theta_1 = theta_1 - M_PI;}

  if (theta_2 < - M_PI_2){theta_2 = theta_2 + M_PI;}
  if (theta_2 > M_PI_2){theta_2 = theta_2 - M_PI;}

        
  if (fabs(theta_1) > fabs(theta_2)){
        	theta_aux = theta_2;
        	theta_2 = theta_1;
        	theta_1 = theta_aux;}

  double fabs_theta_1 = fabs(theta_1);
	beta = 2 * ( fabs_theta_1 + sin ( fabs ( theta_1 + theta_2 )));

  dist = sqrt( pos[1]*pos[1] + pos[0]*pos[0] ) * rescale;

	if ( (dist > 0.0 && dist <= 10.0 && beta < PI_2_69) || (  
	   ( (dist > 0.0 && dist <= 10.0 && beta < PI_1_35) && (fabs_theta_1 < PI_5_9 ) && (fabs(theta_2) < PI_5_9 ) ) ) ){ 
        double dist_1 = 1/dist;
        exp_ = - beta * dist_1 * beta * dist_1 - 2 * powf ( beta * dist_1,7) - dist * dist * 0.011111; // 1/90 = 0.011111
		    return  kappa * powf( euler , exp_);}

	else { return 0.0;}
    }

  nest::TopologyParameter * clone() const
    { return new PlosOne_J(*this); }

private:
  double kappa, orientation_i, orientation_j,rescale;
};




class PlosOne_W: public nest::TopologyParameter
{
public:
  PlosOne_W(const DictionaryDatum& d):
    TopologyParameter(d),
    orientation_i(0.0),
    orientation_j(0.0),
    kappa(1.0),
    rescale(1.0)
    {
      updateValue<double>(d, "orientation_i", orientation_i);
      updateValue<double>(d, "orientation_j", orientation_j);
      updateValue<double>(d, "kappa", kappa);
      updateValue<double>(d, "rescale", rescale);
    }

  double raw_value(const nest::Position<2>& pos,
                   librandom::RngPtr&) const
    { 

  double theta_1,theta_2,theta_aux,beta,dist,exp_1,exp_2,angle; 
  double fabs_or_i_or_j = fabs(orientation_i - orientation_j); 
  double theta_diff = fmin(fabs_or_i_or_j, M_PI - fabs_or_i_or_j);

	angle = atan( pos[1] / pos[0] );

	theta_1 = angle - orientation_i; theta_2 = angle - orientation_j;

  if (theta_1 < - M_PI_2){theta_1 = theta_1 + M_PI;}
  if (theta_1 > M_PI_2){theta_1 = theta_1 - M_PI;}

  if (theta_2 < - M_PI_2){theta_2 = theta_2 + M_PI;}
  if (theta_2 > M_PI_2){theta_2 = theta_2 - M_PI;}


  if (fabs(theta_1) > fabs(theta_2)){
        	theta_aux = theta_2;
        	theta_2 = theta_1;
        	theta_1 = theta_aux;
	}
  double fabs_theta_1 = fabs(theta_1);
  beta = 2*(fabs_theta_1 + sin(fabs(theta_1 + theta_2)));

  dist = sqrt( pos[1] * pos[1] + pos[0] * pos[0] ) * rescale ;

	if ((dist < 0.001 || dist >= 10.0 || beta < PI_1_35) || fabs(theta_diff) >= PI_3 || fabs_theta_1 < PI_11_999){ 
		return 0.0;}
	else {
		exp_1 = -0.4 * pow( beta / dist , 1.5);
 	  exp_2 = -powf( theta_diff * M_4_PI , 1.5);
    return kappa * ( 1 - powf( euler , exp_1) ) * powf( euler , exp_2);
	}
    }

  nest::TopologyParameter * clone() const
    { return new PlosOne_W(*this); }

private:
  double orientation_i, orientation_j, kappa,rescale;
};

} // namespace mynest

#endif 




