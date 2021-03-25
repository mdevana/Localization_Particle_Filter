/**
 * particle_filter.cpp
 *
 * Created on: Dec 12, 2016
 * Author: Tiffany Huang
 */

#include "particle_filter.h"

#include <math.h>
#include <algorithm>
#include <iostream>
#include <iterator>
#include <numeric>
#include <random>
#include <string>
#include <vector>

#include "helper_functions.h"

using std::string;
using std::vector;
using std::normal_distribution;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
  /**
   * TODO: Set the number of particles. Initialize all particles to 
   *   first position (based on estimates of x, y, theta and their uncertainties
   *   from GPS) and all weights to 1. 
   * TODO: Add random Gaussian noise to each particle.
   * NOTE: Consult particle_filter.h for more information about this method 
   *   (and others in this file).
   */
  
  num_particles = 100;  // TODO: Set the number of particles
  std::vector<double> wts(num_particles,1); // initialise a vector of equal weights
  weights = wts;
  
  std::default_random_engine gen;
  
  normal_distribution<double> dist_x(x, std[0]);
  normal_distribution<double> dist_y(y, std[1]);
  normal_distribution<double> dist_theta(theta, std[2]);
  
  std::vector<Particle> P;
  
  for (int i = 0; i < num_particles; ++i) {
	  
	 Particle current_particle;
	 current_particle.id=(i+1);
	 current_particle.x=dist_x(gen);
	 current_particle.y=dist_y(gen);
	 current_particle.theta=dist_theta(gen);
	 current_particle.weight=1;
	 P.push_back(current_particle);
	  
  }
  particles=P;

}

void ParticleFilter::prediction(double delta_t, double std_pos[], 
                                double velocity, double yaw_rate) {
  /**
   * TODO: Add measurements to each particle and add random Gaussian noise.
   * NOTE: When adding noise you may find std::normal_distribution 
   *   and std::default_random_engine useful.
   *  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
   *  http://www.cplusplus.com/reference/random/default_random_engine/
   */
   
   std::default_random_engine gen;
  
   normal_distribution<double> dist_x(0, std_pos[0]);
   normal_distribution<double> dist_y(0, std_pos[1]);
   normal_distribution<double> dist_theta(0, std_pos[2]);
   
   std::cout<< " Velocity = "<< velocity << " Yaw Rate "<< yaw_rate << " delta=  "<<delta_t<<endl;
   
   
   for (int i = 0; i < num_particles; ++i) {
	  
	 //std::cout<<particles[i].id;
	 
	 particles[i].x= ( particles[i].x + (velocity/yaw_rate) * (sin( particles[i].theta + delta_t * yaw_rate) - sin( particles[i].theta)) ) + dist_x(gen);
	 particles[i].y= ( particles[i].y + (velocity/yaw_rate) * (cos( particles[i].theta ) - cos( particles[i].theta +  delta_t * yaw_rate))) + dist_y(gen);
	 particles[i].theta = ( particles[i].theta + ( yaw_rate * delta_t) ) + dist_theta(gen);
	 
	 std::cout<< " X= "<< particles[i].x << " Y= "<< particles[i].y << " theta=  "<<particles[i].theta<<endl;

  }
   
   

}

void ParticleFilter::dataAssociation(vector<LandmarkObs> predicted, 
                                     vector<LandmarkObs>& observations) {
  /**
   * TODO: Find the predicted measurement that is closest to each 
   *   observed measurement and assign the observed measurement to this 
   *   particular landmark.
   * NOTE: this method will NOT be called by the grading code. But you will 
   *   probably find it useful to implement this method and use it as a helper 
   *   during the updateWeights phase.
   */
   double calc_dist;
   for (std::size_t i = 0; i < observations.size(); i++) {
	double min_dist = 10000000;
	LandmarkObs close_LandMarkObs;
	for (std::size_t j = 0; j < predicted.size(); j++){
		calc_dist=dist ( predicted[i].x, predicted[i].y, observations[j].x, observations[j].y);
		if (calc_dist < min_dist) {
			min_dist=calc_dist;
			close_LandMarkObs= predicted[j];
		}
		observations[i].id= close_LandMarkObs.id;
	}
	
   }

}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
                                   const vector<LandmarkObs> &observations, 
                                   const Map &map_landmarks) {
  /**
   * TODO: Update the weights of each particle using a mult-variate Gaussian 
   *   distribution. You can read more about this distribution here: 
   *   https://en.wikipedia.org/wiki/Multivariate_normal_distribution
   * NOTE: The observations are given in the VEHICLE'S coordinate system. 
   *   Your particles are located according to the MAP'S coordinate system. 
   *   You will need to transform between the two systems. Keep in mind that
   *   this transformation requires both rotation AND translation (but no scaling).
   *   The following is a good resource for the theory:
   *   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
   *   and the following is a good resource for the actual equation to implement
   *   (look at equation 3.33) http://planning.cs.uiuc.edu/node99.html
   */
   double head_angle, x_m, y_m,x_c,y_c,mu_x,mu_y;
   
   
   vector<LandmarkObs> Landmarks;
   for (std::size_t l=0; l < map_landmarks.landmark_list.size(); l++) {
	   LandmarkObs new_landmark;
	   new_landmark.id=map_landmarks.landmark_list[l].id_i;
	   new_landmark.x=map_landmarks.landmark_list[l].x_f;
	   new_landmark.y=map_landmarks.landmark_list[l].y_f;
	   Landmarks.push_back(new_landmark);
   }
   
   vector<LandmarkObs> Landmarks_observations=observations; 

   dataAssociation(Landmarks,Landmarks_observations);
   
   for (int i = 0; i < num_particles; ++i) {
	   
	    double final_weight=1.0;
		
		std::cout<< " X= "<< particles[i].x << " Y= "<< particles[i].y << " theta=  "<<particles[i].theta;
		
		for (std::size_t v=0; v < observations.size(); v++) {
			
			head_angle = particles[i].theta; 
			
			x_c = observations[v].x;
			y_c = observations[v].y;
			
			std::cout<< " X_c= "<< x_c << " Y_c = "<< y_c << " theta=  "<<head_angle;
			
			x_m = particles[i].x + ( cos(head_angle) * x_c  - sin(head_angle) * y_c );
			y_m = particles[i].y + ( sin(head_angle) * x_c  - cos(head_angle) * y_c );
			
			for (std::size_t k=0; k < map_landmarks.landmark_list.size(); k++) {
				if ( map_landmarks.landmark_list[k].id_i == observations[v].id ) {
					mu_x = map_landmarks.landmark_list[k].x_f;
					mu_y = map_landmarks.landmark_list[k].y_f;
				}
			}
			
			double d = multi_prob_dist(x_m,y_m,mu_x,mu_y,std_landmark[0],std_landmark[1]);
			std::cout << "multi prob value = " << d << std::endl;
			final_weight*=multi_prob_dist(x_m,y_m,mu_x,mu_y,std_landmark[0],std_landmark[1]);
		}

		
		weights[i] = final_weight;
		particles[i].weight = final_weight;
	}
		
	NormalizeWeights();
	// set normalised Weights to particles
	for (int j = 0; j < num_particles; ++j) {
		particles[j].weight = weights[j];
	}

}

void ParticleFilter::NormalizeWeights() {
	
	double mag=0.0;
	for (std::size_t i=0; i < weights.size(); i++){
		mag+=weights[i];
	}
	double inv_mag= (1.0 / mag);
	for (std::size_t i=0; i < weights.size(); i++){
		weights[i]=weights[i]*inv_mag;
	}
	
}


void ParticleFilter::resample() {
  /**
   * TODO: Resample particles with replacement with probability proportional 
   *   to their weight. 
   * NOTE: You may find std::discrete_distribution helpful here.
   *   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
   */
   
   std::random_device rd;
   std::mt19937 gen(rd());
   
   std::discrete_distribution<> d(0,1);
   
   std::vector<Particle> particles_tmp;
   
   double beta = 0.0;
   double max_weight= *max_element(weights.begin(),weights.end());
   int index = (int) (d(gen) * num_particles);
   
   for( int i = 0;i < num_particles; i++){
	   
	   beta += (d(gen) * 2.0 * max_weight);
	   while ( beta > weights[index] ){
		   beta -= weights[index];
		   index = ( index + 1 ) % num_particles ;
	   }
		particles_tmp.push_back(particles[index]);
	   
   }
   
   particles=particles_tmp;

}

void ParticleFilter::SetAssociations(Particle& particle, 
                                     const vector<int>& associations, 
                                     const vector<double>& sense_x, 
                                     const vector<double>& sense_y) {
  // particle: the particle to which assign each listed association, 
  //   and association's (x,y) world coordinates mapping
  // associations: The landmark id that goes along with each listed association
  // sense_x: the associations x mapping already converted to world coordinates
  // sense_y: the associations y mapping already converted to world coordinates
  particle.associations= associations;
  particle.sense_x = sense_x;
  particle.sense_y = sense_y;
}

string ParticleFilter::getAssociations(Particle best) {
  vector<int> v = best.associations;
  std::stringstream ss;
  copy(v.begin(), v.end(), std::ostream_iterator<int>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length()-1);  // get rid of the trailing space
  return s;
}

string ParticleFilter::getSenseCoord(Particle best, string coord) {
  vector<double> v;

  if (coord == "X") {
    v = best.sense_x;
  } else {
    v = best.sense_y;
  }

  std::stringstream ss;
  copy(v.begin(), v.end(), std::ostream_iterator<float>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length()-1);  // get rid of the trailing space
  return s;
}