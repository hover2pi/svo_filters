from distutils.extension import Extension

def get_package_data():
    return {'svo_filters': ['data/filters/*', 'data/spectra/*']}
