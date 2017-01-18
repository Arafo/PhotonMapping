/*********************************************************************************
Copyright (C) 2014 Adrian Jarabo (ajarabo@unizar.es)
Copyright (C) 2014 Diego Gutierrez (diegog@unizar.es)
All rights reserved.

This is an educational Ray Tracer developed for the course 'Informatica Grafica'
(Computer Graphics) tought at Universidad de Zaragoza (Spain). As such, it does not 
intend to be fast or general, but just to provide an educational tool for undergraduate
students. 

This software is provided as is, and any express or implied warranties are disclaimed.
In no event shall copyright holders be liable for any damage.
**********************************************************************************/

#define _USE_MATH_DEFINES
#include <cmath>

#include "PhotonMapping.h"
#include "World.h"
#include "Intersection.h"
#include "Ray.h"
#include "BSDF.h"

//*********************************************************************
// Compute the photons by tracing the Ray 'r' from the light source
// through the scene, and by storing the intersections with matter
// in the lists 'xx_photons', storing the diffuse (global) and caustic
// photons respectively. For efficiency, both are computed at the same
// time, since computing them separately would result into a lost of
// several samples marked as caustic or diffuse.
// Same goes with the boolean 'direct', that specifies if direct 
// photons (from light to surface) are being stored or not. 
// The initial traced photon has energy defined by the tristimulus
// 'p', that accounts for the emitted power of the light source.
// The function will return true when there are more photons (caustic
// or diffuse) to be shot, and false otherwise.
//---------------------------------------------------------------------
bool PhotonMapping::trace_ray(const Ray& r, const Vector3 &p, 
			   std::list<Photon> &global_photons, std::list<Photon> &caustic_photons, bool direct)
{

	//Check if max number of shots done...
	if( ++m_nb_current_shots > m_max_nb_shots )
	{
		return false;
	}
	
	// Compute irradiance photon's energy
	Vector3 energy(p);
	
	Ray photon_ray(r);
	photon_ray.shift();

	bool is_caustic_particle = false;

	//Iterate the path
	while(1)
	{
		// Throw ray and update current_it
		Intersection it;
		world->first_intersection(photon_ray, it);

		if( !it.did_hit() )
			break;

		//if (photon_ray.get_level() > 0 && direct)
			//break;

		//Check if has hit a delta material...
		if( it.intersected()->material()->is_delta() )
		{
			// If delta material, then is caustic...
			// Don't store the photon!
			is_caustic_particle = true;
		}
		else if (photon_ray.get_level() > 0 || direct)
		{
			//If non-delta material, store the photon!
			if( is_caustic_particle )	
			{				
				//If caustic particle, store in caustics
				if( caustic_photons.size() < m_nb_caustic_photons )
					caustic_photons.push_back( Photon(it.get_position(), photon_ray.get_direction(), energy ));
			}
			else						
			{
				//If non-caustic particle, store in global
				if( global_photons.size() < m_nb_global_photons )
					global_photons.push_back( Photon(it.get_position(), photon_ray.get_direction(), energy ));
			}
			is_caustic_particle = false;
		}	
		
		Real pdf;

		Vector3 surf_albedo = it.intersected()->material()->get_albedo(it);
		Real avg_surf_albedo = surf_albedo.avg();

		Real epsilon2 = static_cast<Real>(rand())/static_cast<Real>(RAND_MAX);
		while (epsilon2 < 0.)
			epsilon2 = static_cast<Real>(rand())/static_cast<Real>(RAND_MAX);
		
		if (epsilon2 > avg_surf_albedo || photon_ray.get_level() > 20 ) 
			break;
			
		// Random walk's next step
		// Get sampled direction plus pdf, and update attenuation
		it.intersected()->material()->get_outgoing_sample_ray(it, photon_ray, pdf );

		// Shade...
		energy = energy*surf_albedo;
		if( !it.intersected()->material()->is_delta() )
			energy *= dot_abs(it.get_normal(), photon_ray.get_direction())/3.14159;		

		energy = energy /(pdf*avg_surf_albedo);
	}
	
	if( caustic_photons.size() == m_nb_caustic_photons && 
		global_photons.size() == m_nb_global_photons )
	{
		m_max_nb_shots = m_nb_current_shots-1;
		return false;
	}

	return true;
}

//*********************************************************************
// TODO: Implement the preprocess step of photon mapping,
// where the photons are traced through the scene. To do it,
// you need to follow these steps for each shoot:
//  1 - Sample a world's light source in the scene to create
//		the initial direct photon from the light source.
//	2 - Trace the photon through the scene storing the inter-
//		sections between the photons and matter. You can use
//		the function 'trace_ray' for this purpose.
//	3 - Finally, once all the photons have been shot, you'll
//		need to build the photon maps, that will be used later
//		for rendering. 
//		NOTE: Careful with function
//---------------------------------------------------------------------
void PhotonMapping::preprocess()
{
	list<Photon> *global_photons = new list<Photon>();
	list<Photon> *caustic_photons = new list<Photon>();

	Ray *rayo;
	Vector3 intensidad;
	Vector3 direccion(0, 0, 0);

	vector<LightSource*> luces = world->light_source_list;
	for (int i = 0; i < luces.size(); i++) {
		LightSource *luz = luces[i];
		intensidad = Vector3(luz->get_intensities() / (m_max_nb_shots /*/ 3.5*/));

		/*
		2 - Trace the photon through the scene storing the inter-
			sections between the photons and matter. You can use
			the function 'trace_ray' for this purpose.
		*/
		do {
			/*
			1 - Sample a world's light source in the scene to create
				the initial direct photon from the light source.
			*/
			Real x, y, z;
			do {
				x = -1 + static_cast <Real> (rand()) / (static_cast <Real> (RAND_MAX / 2));
				y = -1 + static_cast <Real> (rand()) / (static_cast <Real> (RAND_MAX / 2));
				z = -1 + static_cast <Real> (rand()) / (static_cast <Real> (RAND_MAX / 2));
			} while (x * x + y * y + z * z > 1);

			direccion = Vector3(x, y, z);
			rayo = new Ray(luz->get_position(), direccion.normalize(), 0);
			//trace_ray(*rayo, intensidad, *global_photons, *caustic_photons, true);
		} while (trace_ray(*rayo, intensidad, *global_photons, *caustic_photons, false));

		/*
		3 - Finally, once all the photons have been shot, you'll
			need to build the photon maps, that will be used later
			for rendering.
		*/
		for (int j = 0; j < global_photons->size(); j++) {
			Photon photon = global_photons->front();
			global_photons->pop_front();
			m_global_map.store(vector<Real>(photon.position.data, photon.position.data + 3), photon);
		}

		for (int j = 0; j < caustic_photons->size(); j++) {
			Photon photon = caustic_photons->front();
			caustic_photons->pop_front();
			m_caustics_map.store(vector<Real>(photon.position.data, photon.position.data + 3), photon);
		}
	}
	
	if (global_photons->size() > 0)
		m_global_map.balance();
	if (caustic_photons->size() > 0)
		m_caustics_map.balance();
}

//*********************************************************************
// TODO: Implement the function that computes the rendering equation 
// using radiance estimation with photon mapping, using the photon
// maps computed as a proprocess. Note that you will need to handle
// both direct and global illumination, together with recursive the 
// recursive evaluation of delta materials. For an optimal implemen-
// tation you should be able to do it iteratively.
// In principle, the class is prepared to perform radiance estimation
// using k-nearest neighbors ('m_nb_photons') to define the bandwidth
// of the kernel.
//---------------------------------------------------------------------
Vector3 PhotonMapping::shade(Intersection &it0)const
{
	Vector3 L(0);
	Intersection it(it0);

	/*
	* Luz directa
	*/

	Vector3 albedo = it.intersected()->material()->get_albedo(it);
	Vector3 ambiental = world->get_ambient();
	int brillo = 8;

	for (LightSource* luz : world->light_source_list) {
		// Termino ambiental
		L += albedo * ambiental;

		// Intensidad luz
		Vector3 intensidadLuz = luz->get_intensities();

		// Direccion rayo desde luz
		Vector3 direccionLuzAux = luz->get_incoming_direction(it.get_position());	  
		Vector3 direccionLuz = direccionLuzAux * -1;

		// Direccion rayo del punto al ojo
		Vector3 direccionRayoAux = it.get_ray().get_direction();  
		// Direccion rayo del punto al ojo corregido
		Vector3 direccionRayo = direccionRayoAux * -1;

		if (luz->is_visible(it.get_position()) && !it.intersected()->material()->is_delta()) {
			Real lambert = 0;
			
			// Termino difuso
			lambert = it.get_normal().dot(direccionLuz);
			if (lambert > 0) {
				L += lambert * intensidadLuz * albedo;
			}

			// Termino especular
			Real n = it.intersected()->material()->get_specular(it);

			if (n > 0) {
				// 2N(N escalar L) - L
				Vector3 aux = it.get_normal() * 2 * lambert;
				Vector3 LR = aux - direccionLuz;
				Real producto = direccionRayo.dot(LR);
				if (producto > 0) {
 					Real especular = pow(producto, n);
					L += especular  * intensidadLuz * albedo;
				}
			}
		}
	}

	if (it.intersected()->material()->is_delta()) {
		Ray rayo;
		Real ignore(0);
		Intersection itdelta;
		it.intersected()->material()->get_outgoing_sample_ray(it, rayo, ignore);
		rayo.shift();
		world->first_intersection(rayo, itdelta);
		int i = 0;
		while (itdelta.did_hit() && itdelta.intersected()->material()->is_delta() && i < 5) {
			itdelta.intersected()->material()->get_outgoing_sample_ray(itdelta, rayo, ignore);
			world->first_intersection(rayo, itdelta);
			i++;
		}

		if (itdelta.did_hit() && i < 5) {
			Vector3 color = shade(itdelta);
			L += color;
		}
	}

	/*
	* Luz indirecta
	*/

	if (!it.intersected()->material()->is_delta()) {
	vector<const KDTree<Photon, 3>::Node*> global_photons;
	vector<const KDTree<Photon, 3>::Node*> caustic_photons;

	Real maxDistanciaGlobal;
	Real maxDistanciaCaustica;

	// Buscar los fotones globales más cercanos
	m_global_map.find(vector<Real>(it.get_position().data, it.get_position().data + 3), m_nb_photons, global_photons, maxDistanciaGlobal);

	// Buscar los fotones causticos más cercanos
	m_caustics_map.find(vector<Real>(it.get_position().data, it.get_position().data + 3), m_nb_photons, caustic_photons, maxDistanciaCaustica);

	// Cono
	Real k = 1.0f;
	Real areaConoGlobal = 1 / ((1 - (2 / (3 * k))) * M_PI * maxDistanciaGlobal * maxDistanciaGlobal);
	Real areaConoCaustica = 1 /  ((1 - (2 / (3 * k))) * M_PI * maxDistanciaCaustica * maxDistanciaCaustica);

	Vector3 acumuladodGlobal(0, 0, 0);
	Vector3 acumuladodCaustica(0, 0, 0);

	// Estimación de radiancia con el mapa de fotones global
	for (int i = 0; i < global_photons.size(); i++) {
		const KDTree<Photon, 3>::Node* nodo = global_photons.at(i);
		Photon foton = nodo->data();

		Real distancia = (foton.position - it.get_position()).length();
		Real wpc = 1 - (distancia / k * maxDistanciaGlobal);
		acumuladodGlobal += albedo * foton.flux * wpc;
	}

	acumuladodGlobal *= areaConoGlobal;
	L += acumuladodGlobal;
	
	// Estimación de radiancia con el mapa de fotones caustico
	for (int i = 0; i < caustic_photons.size(); i++) {
		const KDTree<Photon, 3>::Node* nodo = caustic_photons.at(i);
		Photon foton = nodo->data();

		Real distancia = (foton.position - it.get_position()).length();
		Real wpc = 1 - (distancia / k * maxDistanciaCaustica);
		acumuladodCaustica += albedo * foton.flux * wpc;
	}

	acumuladodCaustica *= areaConoCaustica;
	L += acumuladodCaustica;
	}

	return L;

	//**********************************************************************
	// The following piece of code is included here for two reasons: first
	// it works as a 'hello world' code to check that everthing compiles 
	// just fine, and second, to illustrate some of the functions that you 
	// will need when doing the work. Goes without saying: remove the 
	// pieces of code that you won't be using.
	//
	/*unsigned int debug_mode = 1;

	switch (debug_mode)
	{
	case 1:
		// ----------------------------------------------------------------
		// Display Albedo Only
		L = it.intersected()->material()->get_albedo(it);
		break;
	case 2:
		// ----------------------------------------------------------------
		// Display Normal Buffer
		L = it.get_normal();
		break;
	case 3:
		// ----------------------------------------------------------------
		// Display whether the material is specular (or refractive) 
		L = Vector3(it.intersected()->material()->is_delta());
		break;

	case 4:
		// ----------------------------------------------------------------
		// Display incoming illumination from light(0)
		L = world->light(0).get_incoming_light(it.get_position());
		break;

	case 5:
		// ----------------------------------------------------------------
		// Display incoming direction from light(0)
		L = world->light(0).get_incoming_direction(it.get_position());
		break;

	case 6:
		// ----------------------------------------------------------------
		// Check Visibility from light(0)
		if (world->light(0).is_visible(it.get_position()))
			L = Vector3(1.);
		break;
	}
	// End of exampled code
	//**********************************************************************

	return L;*/
}