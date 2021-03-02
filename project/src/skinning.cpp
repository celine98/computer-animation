#include "skinning.hpp"

namespace vcl
{
	void normalize_weights(buffer<buffer<float>>& weights)
	{
		size_t const N = weights.size();
		for (size_t k = 0; k < N; ++k) {
			float s = 0.0f;
			for(float w : weights[k]) s += w;
			assert_vcl_no_msg(s>1e-5f);
			for(float& w : weights[k]) w /= s;
		}
	}


	// Linear Blend Skinning
	void skinning_LBS_compute(
		buffer<vec3>& position_skinned,  // position to deform
		buffer<vec3>& normal_skinned,    // normal to deform
		buffer<affine_rt> const& skeleton_current,    // rigid transforms for skeleton joints in current pose
		buffer<affine_rt> const& skeleton_rest_pose,  // rigid transforms of skeleton joints in rest pose
		buffer<vec3> const& position_rest_pose,       // vertex positions of the mesh in rest pose
		buffer<vec3> const& normal_rest_pose,         // normal coordinates of the mesh in rest pose
		rig_structure const& rig)                     // information of the skinning weights (joints and weights associated to a vertex)
	{
		size_t const N_vertex = position_rest_pose.size();
		size_t const N_joint = skeleton_current.size();

		// Sanity check on sizes of buffers
		assert_vcl_no_msg(position_skinned.size()==N_vertex);
		assert_vcl_no_msg(normal_skinned.size()==N_vertex);
		assert_vcl_no_msg(normal_rest_pose.size()==N_vertex);
		assert_vcl_no_msg(skeleton_rest_pose.size()==N_joint);
		assert_vcl_no_msg(rig.joint.size()==N_vertex);
		assert_vcl_no_msg(rig.weight.size()==N_vertex);

		// To do
		//   Compute the Linear Blend Skinning ...

		
		
		for (size_t k = 0; k<N_vertex; k++){
			vec3 const& p0 = position_rest_pose[k];
			vec3 const& n0 = normal_rest_pose[k];
			vec3& p = position_skinned[k];
			vec3& n = normal_skinned[k];
			buffer<int> joints = rig.joint[k];
			buffer<float> weights = rig.weight[k];

			mat4 M_tot;

			for (size_t kj = 0; kj< joints.size(); kj++){
				mat4 T = skeleton_current[joints[kj]].matrix();
				mat4 T0 = inverse(skeleton_rest_pose[joints[kj]]).matrix();
				mat4 Tk = T*T0;

				M_tot += weights[kj]*Tk;

			}

			p = (M_tot*vec4(p0,1.0f)).xyz();
			n = (M_tot*vec4(n0, 0.0f)).xyz();

		}

	}


	void skinning_QBS_compute(
		buffer<vec3>& position_skinned,  // position to deform
		buffer<vec3>& normal_skinned,    // normal to deform
		buffer<affine_rt> const& skeleton_current,    // rigid transforms for skeleton joints in current pose
		buffer<affine_rt> const& skeleton_rest_pose,  // rigid transforms of skeleton joints in rest pose
		buffer<vec3> const& position_rest_pose,       // vertex positions of the mesh in rest pose
		buffer<vec3> const& normal_rest_pose,         // normal coordinates of the mesh in rest pose
		rig_structure const& rig)                     // information of the skinning weights (joints and weights associated to a vertex)
	{

		size_t const N_vertex = position_rest_pose.size();
		size_t const N_joint = skeleton_current.size();

		assert_vcl_no_msg(position_skinned.size()==N_vertex);
		assert_vcl_no_msg(normal_skinned.size()==N_vertex);
		assert_vcl_no_msg(normal_rest_pose.size()==N_vertex);
		assert_vcl_no_msg(skeleton_rest_pose.size()==N_joint);
		assert_vcl_no_msg(rig.joint.size()==N_vertex);
		assert_vcl_no_msg(rig.weight.size()==N_vertex);

		buffer<quaternion> dq;
		dq.resize(N_joint);
		for (size_t kj = 0; kj < N_joint; ++kj) {
			


	}

}

