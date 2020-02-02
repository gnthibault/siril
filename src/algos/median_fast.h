inline __attribute__((always_inline)) float mymin(float a, float b) {
    return b < a ? b : a;
}

inline __attribute__((always_inline))float mymax(float a, float b) {
    return a < b ? b : a;
}

inline __attribute__((always_inline)) float median3x3(float* array)
{
    float tmp = mymin(array[1], array[2]);
    array[2] = mymax(array[1], array[2]);
    array[1] = tmp;
    tmp = mymin(array[4], array[5]);
    array[5] = mymax(array[4], array[5]);
    array[4] = tmp;
    tmp = mymin(array[7], array[8]);
    array[8] = mymax(array[7], array[8]);
    array[7] = tmp;
    tmp = mymin(array[0], array[1]);
    array[1] = mymax(array[0], array[1]);
    array[0] = tmp;
    tmp = mymin(array[3], array[4]);
    array[4] = mymax(array[3], array[4]);
    array[3] = tmp;
    tmp = mymin(array[6], array[7]);
    array[7] = mymax(array[6], array[7]);
    array[6] = tmp;
    tmp = mymin(array[1], array[2]);
    array[2] = mymax(array[1], array[2]);
    array[1] = tmp;
    tmp = mymin(array[4], array[5]);
    array[5] = mymax(array[4], array[5]);
    array[4] = tmp;
    tmp = mymin(array[7], array[8]);
    array[8] = mymax(array[7], array[8]);
    array[3] = mymax(array[0], array[3]);
    array[5] = mymin(array[5], array[8]);
    array[7] = mymax(array[4], tmp);
    tmp = mymin(array[4], tmp);
    array[6] = mymax(array[3], array[6]);
    array[4] = mymax(array[1], tmp);
    array[2] = mymin(array[2], array[5]);
    array[4] = mymin(array[4], array[7]);
    tmp = mymin(array[4], array[2]);
    array[2] = mymax(array[4], array[2]);
    array[4] = mymax(array[6], tmp);
    return mymin(array[4], array[2]);
}
