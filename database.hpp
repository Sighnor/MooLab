#ifndef MOOLAB_DATABASE
#define MOOLAB_DATABASE

#include "controller.hpp"
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
    Character_Controller &cc, 
    int contact, 
    vec3 ltoe, 
    vec3 rtoe, 
    int N = 50, 
    int dN = 10, 
    float weight_vel = 1.f, 
    float weight_avel = 10.f, 
    float weight_toe = 3.f)
{
    int best_frame = 0;
    float best_delta = 10000.f;
    for (int i = 0; i < db.nframes() - 1.5 * N; i++)
    {
        if(contact == db.contacts(i))
        {
            vec3 db_vel = (db.vels(i, 1) + db.vels(i, 2) + db.vels(i, 3) + db.vels(i, 4) + db.vels(i, 5)) / 5.f;
            vec3 db_avel = (db.avels(i, 1) + db.avels(i, 2) + db.avels(i, 3) + db.avels(i, 4) + db.avels(i, 5)) / 5.f;

            float delta_vel0 = length(cc.vel0) < 0.1f && length(db.vels(i, 0)) < 0.1f? 0.f : length(cc.vel0 - db.vels(i, 0)) / length(cc.vel0);
            float delta_vel1 = length(cc.vel1) < 0.1f && length(db.vels(i, 1)) < 0.1f? 0.f : length(cc.vel1 - db.vels(i, 1)) / length(cc.vel1);
            float delta_vel2 = length(cc.vel2) < 0.1f && length(db.vels(i, 2)) < 0.1f? 0.f : length(cc.vel2 - db.vels(i, 2)) / length(cc.vel2);
            float delta_vel3 = length(cc.vel3) < 0.1f && length(db.vels(i, 3)) < 0.1f? 0.f : length(cc.vel3 - db.vels(i, 3)) / length(cc.vel3);
            float delta_vel4 = length(cc.vel4) < 0.1f && length(db.vels(i, 4)) < 0.1f? 0.f : length(cc.vel4 - db.vels(i, 4)) / length(cc.vel4);
            float delta_vel5 = length(cc.vel5) < 0.1f && length(db.vels(i, 5)) < 0.1f? 0.f : length(cc.vel5 - db.vels(i, 5)) / length(cc.vel5);

            float delta_avel0 = abs(cc.avel0.y - db.avels(i, 0).y);
            float delta_avel1 = abs(cc.avel1.y - db.avels(i, 1).y);
            float delta_avel2 = abs(cc.avel2.y - db.avels(i, 2).y);
            float delta_avel3 = abs(cc.avel3.y - db.avels(i, 3).y);
            float delta_avel4 = abs(cc.avel4.y - db.avels(i, 4).y);
            float delta_avel5 = abs(cc.avel5.y - db.avels(i, 5).y);

            float delta_ltoe = length(ltoe - db.features(i, 0));
            float delta_rtoe = length(rtoe - db.features(i, 1));

            float delta_vel = delta_vel0 + delta_vel1 + delta_vel2 + delta_vel3 + delta_vel4 + delta_vel5;
            float delta_avel = delta_avel0 + delta_avel1 + delta_avel2 + delta_avel3 + delta_avel4 + delta_avel5;
            float delta_toe = delta_ltoe + delta_rtoe;

            float temp_delta = weight_vel * delta_vel + weight_avel * delta_avel + weight_toe * delta_toe;

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