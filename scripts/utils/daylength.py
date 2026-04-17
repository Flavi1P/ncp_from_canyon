import numpy as np
from datetime import datetime

def daylength(year, month, day, lat):
    """
    Computes the length of the day (the time between sunrise and sunset) 
    given the year, month, day, and latitude of the location.
    Function uses the Brock model for the computations.

    For more information see, for example:
    Forsythe et al., "A model comparison for daylength as a function of 
    latitude and day of year", Ecological Modelling, 1995.

    Parameters
    ----------
    year : int
        Year (e.g., 2024)
    month : int
        Month of the year (1 for January, 12 for December)
    day : int
        Day of the month (1 to 31, depending on the month)
    lat : float
        Latitude of the location in degrees. Positive values
        for north and negative for south.

    Returns
    -------
    d : float
        Daylength in hours.
    """
    # Compute day of the year from the given date
    date = datetime(year, month, day)
    dayOfYear = date.timetuple().tm_yday

    # Convert latitude to radians
    latInRad = np.deg2rad(lat)

    # Compute declination of the Earth
    declinationOfEarth = 23.45 * np.sin(np.deg2rad(360.0 * (283.0 + dayOfYear) / 365.0))

    # Handle edge cases where there is constant daylight or darkness
    if -np.tan(latInRad) * np.tan(np.deg2rad(declinationOfEarth)) <= -1.0:
        return 24.0
    elif -np.tan(latInRad) * np.tan(np.deg2rad(declinationOfEarth)) >= 1.0:
        return 0.0
    else:
        # Calculate hour angle
        hourAngle = np.rad2deg(np.arccos(-np.tan(latInRad) * np.tan(np.deg2rad(declinationOfEarth))))
        return 2.0 * hourAngle / 15.0
