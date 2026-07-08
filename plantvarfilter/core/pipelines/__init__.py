"""
PlantOmicsGWAS Core Pipelines

Intentionally empty of imports. Each pipeline step module is imported
lazily and individually by plantvarfilter.core.pipeline_factory,
only when that specific step is actually used. Importing this
package must stay cheap and dependency-free.
"""