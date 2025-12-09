## Cleaned Rat Dataset

**`Unique Key`**

- `[int64]` Unique identifier of rat sighting.

**`Incident Zip`**

- `[int64]` Zip code of rat sighting.


## Cleaned Restaurant Dataset

**`ZIPCODE`**

- `[int64]` Zip code of restaurant.

**`VIOLATION DESCRIPTION`**

- `[object]` Description of violations found during inspection.

**`CRITICAL FLAG`**

- `[object]` Whether restaurant inspection condition was determined to be critical or not.

**`SCORE`**

- `[float64]` Scored between 1-100 on restaurant quality, where lower scores mean less violations (better),


## Merged Dataset

**`Incident Zip`**

- `[int64]` Individual NYC zip codes found in both datasets.

**`ZIPCODE`**

- `[float64]` Individual NYC zip codes found in both datasets.

**`rat_sightings`**

- `[int64]` Total number of rat sightings counted per zip code.

**`total_inspections`**

- `[float64]` Total number of restaurant inspections counted per zip code.

**`unique_restaurants`**

- `[float64]` Total number of unique restaurants counted in the NYC Restaurant Inspections dataset.

**`avg_score`**

- `[float64]` Average score given to restaurants per zip code.

**`median_score`**

- `[float64]` Median score given to restaurants per zip code.

**`critical_rate`**

- `[float64]` Average number of restaurants that were flagged as "Critical" during restaurant inspection, per zip code.
