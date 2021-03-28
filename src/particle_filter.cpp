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

static std::default_random_engine gen;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
  /**
   * TODO: Set the number of particles. Initialize all particles to 
   *   first position (based on estimates of x, y, theta and their uncertainties
   *   from GPS) and all weights to 1. 
   * TODO: Add random Gaussian noise to each particle.
   * NOTE: Consult particle_filter.h for more information about this method 
   *   (and others in this file).
   */
  
  if (is_initialized==true){
	  return;
  }
  
  num_particles = 1;  // TODO: Set the number of particles
  std::vector<double> wts(num_particles,1); // initialise a vector of equal weights
  weights = wts;
  
  
  
  normal_distribution<double> dist_x(x, std[0]);
  normal_distribution<double> dist_y(y, std[1]);
  normal_distribution<double> dist_theta(theta, std[2]);
  
  //std::vector<Particle> P;
  
  //std::cout<< " initial X = "<< x << " initial Y "<< y <<std::endl;
  //std::cout<< " dist X = "<< dist_x(gen) << " dist Y "<< dist_y(gen) <<std::endl;
  
  
  for (int i = 0; i < num_particles; ++i) {
	  
	 Particle current_particle;
	 current_particle.id=(i+1);
	 current_particle.x = dist_x(gen);
	 current_particle.y = dist_y(gen);
	 //current_particle.x = x;
	 //current_particle.y = y;
	 current_particle.theta=  dist_theta(gen);
	 //current_particle.theta=  theta;
	 current_particle.weight=1.0;
	 particles.push_back(current_particle);
	  
  }
  
  is_initialized=true;

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
   
   
   
   normal_distribution<double> dist_x(0, std_pos[0]);
   normal_distribution<double> dist_y(0, std_pos[1]);
   normal_distribution<double> dist_theta(0, std_pos[2]);
  
   
   //std::cout<<" Prediction step"<<std::endl;
   
   //std::cout<< " Velocity = "<< velocity << " Yaw Rate "<< yaw_rate << " delta=  "<<delta_t<<std::endl;
   
   
   for (int i = 0; i < num_particles; ++i) {
	  
	 //std::cout<<particles[i].id;
	 
	 
	 //std::cout<<"before prediction" <<" X= "<< particles[i].x << " Y= "<< particles[i].y << " theta=  "<<particles[i].theta<<std::endl;
	 
	 if (fabs(yaw_rate)< 0.00001){ // no rotation motion
	 
		 particles[i].x = ( particles[i].x + velocity *  delta_t * cos(particles[i].theta) ) + dist_x(gen);
		 particles[i].y = ( particles[i].y + velocity *  delta_t * sin(particles[i].theta) ) + dist_y(gen);
		 particles[i].theta = particles[i].theta + dist_theta(gen);
		 
	 }
	 else {
		 
		 particles[i].x = ( particles[i].x + (velocity/yaw_rate) * (sin( particles[i].theta + delta_t * yaw_rate) - sin( particles[i].theta)) ) + dist_x(gen);
		 particles[i].y = ( particles[i].y + (velocity/yaw_rate) * (cos( particles[i].theta ) - cos( particles[i].theta +  delta_t * yaw_rate))) + dist_y(gen);
		 particles[i].theta = ( particles[i].theta + ( yaw_rate * delta_t) ) + dist_theta(gen);
		 
		 
		 
		 //particles[i].x = ( particles[i].x + (velocity/yaw_rate) * (sin( particles[i].theta + delta_t * yaw_rate) - sin( particles[i].theta)) ) ;
		 //particles[i].y = ( particles[i].y + (velocity/yaw_rate) * (cos( particles[i].theta ) - cos( particles[i].theta +  delta_t * yaw_rate))) ;
		 //particles[i].theta = ( particles[i].theta + ( yaw_rate * delta_t) ) ;
	 }
	 
	 
	 
	 //std::cout<< " X= "<< particles[i].x << " Y= "<< particles[i].y << " theta=  "<<particles[i].theta<<std::endl;

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
	double min_dist = std::numeric_limits<double>::max();
	//std::cout<< " observation X= "<< observations[i].x<< " observation Y= "<< observations[i].y<<std::endl;
	int close_LandMarkObs_id;
	for (std::size_t j = 0; j < predicted.size(); j++){
		calc_dist=dist ( predicted[j].x, predicted[j].y, observations[i].x, observations[i].y);
		//std::cout<< " landmark X= "<< predicted[j].x << " landmark Y= "<< predicted[j].y << " landmark id =  "<<predicted[j].id<< " calc dist ="<<calc_dist<<std::endl;
		if (calc_dist < min_dist) {
			//std::cout<< "min dist "<< min_dist<< " calc dist "<< calc_dist<<std::endl;
			min_dist=calc_dist;
			close_LandMarkObs_id= predicted[j].id;
		}
	}
	//std::cout<<"assigned landmark id :"<<close_LandMarkObs_id<<std::endl;	
	observations[i].id= close_LandMarkObs_id;
		
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
   double head_angle, x_p, y_p, x_m, y_m, x_c,y_c,mu_x,mu_y;
   double calc_dist;
   
   
   for (int i = 0; i < num_particles; ++i) {
	   
	    // Collect data for transformation
		head_angle = particles[i].theta; 
		x_p = particles[i].x;
		y_p = particles[i].y;
		//std::cout<< " X= "<< particles[i].x << " Y= "<< particles[i].y << " theta=  "<<head_angle;
		
		// include only Landmarks inside Sensor Range
		vector<LandmarkObs> Landmarks;
		
		 

		for (std::size_t l=0; l < map_landmarks.landmark_list.size(); l++) {
			LandmarkObs new_landmark;
			new_landmark.id=map_landmarks.landmark_list[l].id_i;
			new_landmark.x=map_landmarks.landmark_list[l].x_f;
			new_landmark.y=map_landmarks.landmark_list[l].y_f;
			//std::cout<< " landmark X= "<< map_landmarks.landmark_list[l].x_f << " landmark Y= "<< map_landmarks.landmark_list[l].y_f << " landmark id =  "<<map_landmarks.landmark_list[l].id_i;
			calc_dist=dist ( particles[i].x , particles[i].y, new_landmark.x, new_landmark.y);
			if (calc_dist <= sensor_range)
				Landmarks.push_back(new_landmark);
		}
		
		vector<LandmarkObs> Landmarks_observations;
		// Loop to make transformation from Car to Map coordinate using particle heading and Position
		for (std::size_t v=0; v < observations.size(); v++) {

			x_c = observations[v].x;
			y_c = observations[v].y;
			
			//std::cout<< " X_c= "<< x_c << " Y_c = "<< y_c << " theta=  "<<head_angle<<std::endl;
			
			x_m = x_p + ( cos(head_angle) * x_c  - sin(head_angle) * y_c );
			y_m = y_p + ( sin(head_angle) * x_c  + cos(head_angle) * y_c );
			
			//std::cout<< " X_p= "<< x_p << " Y_p = "<< y_p << " theta=  "<<head_angle<<std::endl;
			
			//std::cout<< " X_m= "<< x_m << " Y_m = "<< y_m << " cos(headangle) =  "<<cos(head_angle)<<" sin(headangle) =  "<<sin(head_angle)<<std::endl;
			
			Landmarks_observations.push_back(LandmarkObs{observations[v].id, x_m, y_m});
		}
		
		// Associate Observed Landmarks to Landmarks
		dataAssociation(Landmarks,Landmarks_observations);
		
		/*for (std::size_t t=0; t < Landmarks_observations.size(); t++) {
			std::cout<< " landmark_obs X= "<< Landmarks_observations[t].x << " landmark_obs Y= "<< Landmarks_observations[t].y << " landmarkobs id =  "<<Landmarks_observations[t].id;
	   
		}*/
		
		// Calculate Weights using multivariante Gaussian Distribution
		
		
		double final_weight=1.0;
		
		for (std::size_t v=0; v < Landmarks_observations.size(); v++) {
			
				x_m = Landmarks_observations[v].x;
				y_m = Landmarks_observations[v].y;
			
				for (std::size_t k=0; k < Landmarks.size(); k++) {

						if ( Landmarks[k].id == Landmarks_observations[v].id ) {
							mu_x = Landmarks[k].x;
							mu_y = Landmarks[k].y;
							break;
						}
					
				}

				//std::cout<< " X_m= "<< x_m << " Y_m= "<< y_m << " Mu_x=  "<<mu_x<< " Mu_y=  "<<mu_y<< " std_x=  "<<std_landmark[0]<<" std_y=  "<<std_landmark[1]<<std::endl;

				double wt = multi_prob_dist(x_m,y_m,mu_x,mu_y,std_landmark[0],std_landmark[1]);
				if (wt == 0 ){
					wt=0.00001;
				}
				
				//std::cout << "multi prob value = " << wt << std::endl;
				final_weight*=wt;
			
		}

		//std::cout<<" Final calculated weight for particle :"<< final_weight;
		
		particles[i].weight = final_weight;
	}
		
	

}

void ParticleFilter::NormalizeWeights() {
	
	/*double mag=0.0;
	for (std::size_t i=0; i < weights.size(); i++){
		mag+=weights[i];
	}
	double inv_mag= (1.0 / mag);
	for (std::size_t i=0; i < weights.size(); i++){
		weights[i]=weights[i]*inv_mag;
	}*/
	const double A_PI = 3.14159265358979323846;
	for (int i = 0; i < num_particles; ++i) {
		std::cout<<"Particle x Coordinate :" << particles[i].x<<std::endl;
		std::cout<<"Particle y Coordinate :"<< particles[i].y<<std::endl;
		std::cout<<"Particle theta Coordinate :"<< particles[i].theta<<std::endl;
	}
	double delta_t=0.1;
	double std_pos[3]={0.3, 0.3, 0.01};
	double velocity = 110.0; 
	double yaw_rate = A_PI / 8;
	
	//prediction(delta_t,std_pos,velocity,yaw_rate);
	
	for (int i = 0; i < num_particles; ++i) {
		std::cout<<"Particle x Coordinate :" << particles[i].x;
		std::cout<<"Particle y Coordinate :"<< particles[i].y;
		std::cout<<"Particle theta Coordinate :" << particles[i].theta;
	}
	
	
	LandmarkObs OBS1,OBS2,OBS3;
	
	OBS1.x = 2;
	OBS1.y = 2;
	
	OBS2.x = 3;
	OBS2.y = -2;
	
	OBS3.x = 0;
	OBS3.y = -4;
	
	vector<LandmarkObs> obs{OBS1,OBS2,OBS3};
	
	Map map;

    // Declare single_landmark
    Map::single_landmark_s single_landmark_temp;

    // Set values
    single_landmark_temp.id_i = 1;
    single_landmark_temp.x_f  = 5;
    single_landmark_temp.y_f  = 3;

    // Add to landmark list of map
    map.landmark_list.push_back(single_landmark_temp);
	
	single_landmark_temp.id_i = 2;
    single_landmark_temp.x_f  = 2;
    single_landmark_temp.y_f  = 1;

    // Add to landmark list of map
    map.landmark_list.push_back(single_landmark_temp);
	
	single_landmark_temp.id_i = 3;
    single_landmark_temp.x_f  = 6;
    single_landmark_temp.y_f  = 1;

    // Add to landmark list of map
    map.landmark_list.push_back(single_landmark_temp);
	
	single_landmark_temp.id_i = 4;
    single_landmark_temp.x_f  = 7;
    single_landmark_temp.y_f  = 4;

    // Add to landmark list of map
    map.landmark_list.push_back(single_landmark_temp);
	
	single_landmark_temp.id_i = 5;
    single_landmark_temp.x_f  = 4;
    single_landmark_temp.y_f  = 7;

    // Add to landmark list of map
    map.landmark_list.push_back(single_landmark_temp);
	
	double sensor_range = 50;
	double sigma_landmark [2] = {0.3, 0.3};
	
	updateWeights(sensor_range,sigma_landmark,obs,map);
	resample();
	
		
	
}


void ParticleFilter::resample() {
  /**
   * TODO: Resample particles with replacement with probability proportional 
   *   to their weight. 
   * NOTE: You may find std::discrete_distribution helpful here.
   *   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
   */
   
   double max_weight = std::numeric_limits<double>::min();
   
   std::cout<<" in resampling code"<<std::endl;
   
   for( int i = 0;i < num_particles; i++){
	   weights[i]=particles[i].weight;
	   if (particles[i].weight > max_weight)
		max_weight = particles[i].weight;
	   
   }
   
   std::uniform_real_distribution<double> distForBeta(0.0, 2 * max_weight);
   std::uniform_int_distribution<int> distForIndex(0, num_particles - 1);
   
   std::vector<Particle> particles_resampled;
   
   double beta = 0.0;
   
   int index = distForIndex(gen);
   std::cout<<" in resampling code index "<<index <<std::endl;
   
   for( int i = 0; i < num_particles; i++){
	   
	   beta += distForBeta(gen);
	   while ( beta > weights[index] ){
		   beta -= weights[index];
		   index = ( index + 1 ) % num_particles ;
	   }
		particles_resampled.push_back(particles[index]);
	   
   }
   
   particles=particles_resampled;

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