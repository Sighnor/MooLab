#ifndef ENGINE_DATABASE
#define ENGINE_DATABASE

#include "motion.hpp"

struct Database
{
    array1d<int> frames;
    array1d<float> velocities;
    // vel avel ltoe rtoe
    array2d<vec3> features;

    int nframes() const { return frames.size; }
    int nfeatures() const { return features.cols; }
};

void Database_load(Database &db, const char* file_name)
{
    FILE* f = fopen(file_name, "rb");
    assert(f != NULL);

    array1d_read(db.frames, f);
    array1d_read(db.velocities, f);
    array2d_read(db.features, f);

    fclose(f);
}

void Database_save(Database &db, const char* file_name)
{
    FILE* f = fopen(file_name, "wb");
    assert(f != NULL);

    array1d_write(db.frames, f);
    array1d_write(db.velocities, f);
    array2d_write(db.features, f);

    fclose(f);
}

void Database_sort(Database &db, array2d<vec3> features)
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

int Database_search(
    Database &db, 
    vec3 vel0, 
    vec3 avel0,
    vec3 ltoe, 
    vec3 rtoe)
{
    int best_frame = 1000;
    float best_delta_toe = 1.f;
    for (int i = 0; i < db.nframes() - 100; i++)
    {
        float delta_vel0_length = abs(length(vel0) / length(db.features(i, 0)) - 1.f);
        float delta_vel0_dir = dot(vel0, db.features(i, 0)) / length(vel0) / length(db.features(i, 0));
        float delta_vel1_length = abs(length(vel0) / length(db.features(i + 20, 0)) - 1.f);
        float delta_vel1_dir = dot(vel0, db.features(i + 20, 0)) / length(vel0) / length(db.features(i + 20, 0));
        float delta_vel2_length = abs(length(vel0) / length(db.features(i + 40, 0)) - 1.f);
        float delta_vel2_dir = dot(vel0, db.features(i + 40, 0)) / length(vel0) / length(db.features(i + 40, 0));
        float delta_vel3_length = abs(length(vel0) / length(db.features(i + 60, 0)) - 1.f);
        float delta_vel3_dir = dot(vel0, db.features(i + 60, 0)) / length(vel0) / length(db.features(i + 60, 0));
        float delta_vel4_length = abs(length(vel0) / length(db.features(i + 80, 0)) - 1.f);
        float delta_vel4_dir = dot(vel0, db.features(i + 80, 0)) / length(vel0) / length(db.features(i + 80, 0));

        float delta_avel0 = abs(avel0.y - db.features(i, 1).y);
        float delta_avel1 = abs(avel0.y - db.features(i + 20, 1).y);
        float delta_avel2 = abs(avel0.y - db.features(i + 40, 1).y);
        float delta_avel3 = abs(avel0.y - db.features(i + 60, 1).y);
        float delta_avel4 = abs(avel0.y - db.features(i + 80, 1).y);

        float delta_ltoe = length(ltoe - db.features(i, 2));
        float delta_rtoe = length(rtoe - db.features(i, 3));

        if(delta_vel0_dir > 0.9f && 
           delta_vel0_length < 0.1f &&
           delta_avel0 < 1.f; 
           delta_vel1_dir > 0.9f && 
           delta_vel2_dir > 0.9f && 
           delta_vel3_dir > 0.9f && 
           delta_vel4_dir > 0.9f)
        {
            float delta_toe = delta_ltoe + delta_rtoe;
            if(delta_toe < best_delta_toe)
            {
                best_delta_toe = delta_toe;
                best_frame = db.frames(i);
            }
        }
    }
    std::cout << best_delta_toe << ',' << best_frame << std::endl;
    return best_frame;
}

#endif