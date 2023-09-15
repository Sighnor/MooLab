#ifndef ENGINE_DATABASE
#define ENGINE_DATABASE

#include "motion.hpp"

struct Database
{
    array1d<int> frames;
    array1d<int> contacts;
    // frame 0.0 * N 0.2 * N 0.4 * N 0.6 * N 0.8 * N 1.0 * N
    array2d<vec3> vels;
    array2d<vec3> avels;
    // ltoe rtoe
    array2d<vec3> features;

    int nframes() const { return frames.size; }
    int nfeatures() const { return features.cols; }
};

void Database_load(Database &db, const char* file_name)
{
    FILE* f = fopen(file_name, "rb");
    assert(f != NULL);

    array1d_read(db.frames, f);
    array1d_read(db.contacts, f);
    array2d_read(db.vels, f);
    array2d_read(db.avels, f);
    array2d_read(db.features, f);

    fclose(f);
}

void Database_save(Database &db, const char* file_name)
{
    FILE* f = fopen(file_name, "wb");
    assert(f != NULL);

    array1d_write(db.frames, f);
    array1d_write(db.contacts, f);
    array2d_write(db.vels, f);
    array2d_write(db.avels, f);
    array2d_write(db.features, f);

    fclose(f);
}

void Database_sort(Database &db, const slice2d<vec3> features)
{
    db.features.resize(features.rows, features.cols);
    for(int i = 0; i < db.nframes(); i++)
    {
        int frame_id = db.frames(i);
        for(int j = 0; j < db.nfeatures(); j++)
        {
            db.features(i, j) = features(frame_id, j);
        }
    }
}

void Database_revise(Database &db, int dN = 10)
{
    for(int i = 0; i < db.nframes(); i++)
    {
        for(int j = 5; j > 1; j--)
        {
            db.vels(i, j) = vec_to_vel(vel_to_vec(db.vels(i, j - 1), dN * (j - 1) * 0.0166667f), vel_to_vec(db.vels(i, j), dN * j * 0.0166667f), dN * 0.0166667f);
            db.avels(i, j) = quat_to_avel(avel_to_quat(db.avels(i, j - 1), dN * (j - 1) * 0.0166667f), avel_to_quat(db.avels(i, j), dN * j * 0.0166667f), dN * 0.0166667f);
        }
    }
}

int Database_search(
    Database &db, 
    int contact, 
    vec3 vel0, 
    vec3 avel0,
    vec3 vel5, 
    vec3 avel5, 
    vec3 ltoe, 
    vec3 rtoe, 
    int N = 50, 
    int dN = 10)
{
    // vec3 vel1 = vel0 + 0.2f * vel5;
    // vec3 vel2 = vel0 + 0.4f * vel5;
    // vec3 vel3 = vel0 + 0.6f * vel5;
    // vec3 vel4 = vel0 + 0.8f * vel5;

    vec3 vel1 = vel5;
    vec3 vel2 = vel5;
    vec3 vel3 = vel5;
    vec3 vel4 = vel5;

    vec3 avel1 = avel5;
    vec3 avel2 = avel5;
    vec3 avel3 = avel5;
    vec3 avel4 = avel5;

    int best_frame = 0;
    float best_delta = 10000.f;
    for (int i = 0; i < db.nframes() - 1.5 * N; i++)
    {
        float delta_vel0_length = abs(length(vel0) / length(db.vels(i, 0)) - 1.f);
        float delta_vel0_dir = dot(vel0, db.vels(i, 0)) / length(vel0) / length(db.vels(i, 0));
        float delta_vel1_length = abs(length(vel1) / length(db.vels(i, 1)) - 1.f);
        float delta_vel1_dir = dot(vel1, db.vels(i, 1)) / length(vel1) / length(db.vels(i, 1));
        float delta_vel2_length = abs(length(vel2) / length(db.vels(i, 2)) - 1.f);
        float delta_vel2_dir = dot(vel2, db.vels(i, 2)) / length(vel2) / length(db.vels(i, 2));
        float delta_vel3_length = abs(length(vel3) / length(db.vels(i, 3)) - 1.f);
        float delta_vel3_dir = dot(vel3, db.vels(i, 3)) / length(vel3) / length(db.vels(i, 3));
        float delta_vel4_length = abs(length(vel4) / length(db.vels(i, 4)) - 1.f);
        float delta_vel4_dir = dot(vel4, db.vels(i, 4)) / length(vel4) / length(db.vels(i, 4));
        float delta_vel5_length = abs(length(vel5) / length(db.vels(i, 5)) - 1.f);
        float delta_vel5_dir = dot(vel5, db.vels(i, 5)) / length(vel5) / length(db.vels(i, 5));

        float delta_avel0 = abs(avel0.y - db.avels(i, 0).y);
        float delta_avel1 = abs(avel1.y - db.avels(i, 1).y);
        float delta_avel2 = abs(avel2.y - db.avels(i, 2).y);
        float delta_avel3 = abs(avel3.y - db.avels(i, 3).y);
        float delta_avel4 = abs(avel4.y - db.avels(i, 4).y);
        float delta_avel5 = abs(avel5.y - db.avels(i, 5).y);

        float delta_ltoe = length(ltoe - db.features(i, 0));
        float delta_rtoe = length(rtoe - db.features(i, 1));

        float delta_vel_dir = (1.f - delta_vel1_dir) + (1.f - delta_vel3_dir) + (1.f - delta_vel5_dir);
        float delta_avel = delta_avel1 + delta_avel3 + delta_avel5;

        if(contact == db.contacts(i))
        {
            float delta_vel_length = delta_vel1_length + delta_vel3_length + delta_vel5_length;// + delta_vel4_length + delta_vel5_length;
            float delta_toe = delta_ltoe + delta_rtoe;

            float temp_delta = 5.f * delta_vel_dir + 1.f * delta_vel_length + 3.f * delta_toe + 10.f * delta_avel;
            if(temp_delta < best_delta)
            {
                best_delta = temp_delta;
                best_frame = db.frames(i);
            }
        }
    }
    std::cout << best_delta << ',' << best_frame << std::endl;
    return best_frame;
}

#endif