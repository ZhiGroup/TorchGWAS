"""Comparable timing boundaries for streamed calculator validation."""

SETUP_SCAN_WRITE_METRIC = 'setup_scan_and_write_seconds'
SETUP_SCAN_WRITE_BOUNDARY = (
    'After input QC through writer return, including shared preprocessing, '
    'per-device or per-tile setup, scan, iterator cleanup and requested output '
    'publication; excludes input opening/QC and subsequent run metadata. '
    'Publication is durable only when output fsync is enabled.'
)


def streaming_timing(prep_done, write_started, write_finished):
    """Keep the writer-only interval and add a setup-inclusive executor clock.

    A single-device reduction prepares its factor before writer entry. Workers
    in a multi-device reduction prepare factors after writer entry. Comparing
    writer-only times would therefore charge setup only to the latter.
    All fields use one captured writer-return endpoint; these are overlapping
    diagnostics, not additional additive entries in ``phase_seconds``.
    """
    if not prep_done <= write_started <= write_finished:
        raise ValueError('Streaming timing boundaries must be ordered')
    return {
        'scan_setup_seconds': write_started - prep_done,
        'scan_and_write_seconds': write_finished - write_started,
        SETUP_SCAN_WRITE_METRIC: write_finished - prep_done,
        'setup_scan_and_write_boundary': SETUP_SCAN_WRITE_BOUNDARY,
    }
