"""Frozen pre-change trait-coordinate rebase for paired correctness controls."""
import numpy as np
ORIGINAL_FUNCTION_SHA256 = 'e46cf5c54f97beb2c8420c6512cd2fa909eae475cbee1ca1d7e2bbe3912e296b'

def legacy_trait_blocked_significant_chunks(scan_once, significance, n_traits,
                                      trait_block, df, devices=None, queue_depth=None,
                                      partition_context=None):
    """Significant (variant, trait) pairs for all K traits, block by block.

    `reduce='significant'` cannot use `_trait_blocked_reduced_chunks`: that
    one merges a `VariantReduction`'s one-row-per-variant result and expects
    a sixth tuple element carrying the winning trait index. Significance has
    no winner and no merge -- it emits every pair over a threshold, so the
    number of rows per chunk varies and blocks simply concatenate. Routing
    it through the reduced driver raised `not enough values to unpack
    (expected 6, got 5)` the moment blocking first engaged.

    **The threshold is Bonferroni over the WHOLE trait axis, never the
    block.** A block is an implementation detail chosen from the available
    memory, and if it changed which pairs were significant then the results
    would depend on the size of the card they ran on. `n_traits` here is
    always the full K.

    Trait indices come back numbered within their block and are rebased onto
    the full axis before they are yielded, which is the only bookkeeping the
    concatenation needs.
    """
    from torchgwas.linear import _significant_pairs_iterator

    if trait_block < 1:
        raise ValueError("trait_block must be positive")
    if partition_context is not None and not callable(partition_context):
        raise ValueError('Callable indexed producer binding required')
    blocks = [(offset, min(trait_block, n_traits - offset))
              for offset in range(0, n_traits, trait_block)]

    def emit(offset, width, device):
        partition=None if partition_context is None else partition_context(offset,width,device)
        source = scan_once(offset, width, device)
        pairs = _significant_pairs_iterator(source, significance, n_traits, df)
        try:
            for start, end, variant_index, trait_index, beta, t_stat, row_df in pairs:
                item=(start,end,variant_index,trait_index.astype(np.int64)+offset,beta,t_stat,row_df)
                if partition is not None:
                    from torchgwas.sumstats_indexed import PartitionedIndexedChunk
                    item=PartitionedIndexedChunk(item,partition)
                yield item
        finally:
            try:
                if hasattr(pairs, "close"):
                    pairs.close()
            finally:
                if hasattr(source, "close"):
                    source.close()
    if not devices or len(devices) <= 1:
        for offset, width in blocks:
            yield from emit(offset, width, None)
        return

    # One thread per device, each taking its own blocks. Rows are handed
    # back through a bounded queue rather than accumulated: the whole point
    # of this mode is that the output is small but unbounded in principle,
    # so it streams to the writer instead of being held.
    import queue
    import threading

    if queue_depth is None:
        queue_depth = 4 * len(devices)
    if isinstance(queue_depth, bool) or not isinstance(queue_depth, int) or queue_depth < 1:
        raise ValueError('queue_depth must be a positive integer')
    results: queue.Queue = queue.Queue(maxsize=queue_depth)
    failures: list[BaseException] = []
    done = object()
    stop = threading.Event()

    def publish(item):
        while not stop.is_set():
            try:
                results.put(item, timeout=0.05)
                return True
            except queue.Full:
                pass
        return False

    def run(device, assigned):
        iterator = None
        try:
            for offset, width in assigned:
                if stop.is_set():
                    break
                iterator = emit(offset, width, device)
                for item in iterator:
                    if not publish(item):
                        return
                iterator.close()
                iterator = None
        except BaseException as exc:
            failures.append(exc)
            stop.set()
        finally:
            try:
                if iterator is not None:
                    iterator.close()
            except BaseException as exc:
                failures.append(exc)
                stop.set()
            publish(done)

    workers = []
    for index, device in enumerate(devices):
        assigned = blocks[index::len(devices)]
        if not assigned:
            continue
        thread = threading.Thread(target=run, args=(device, assigned),
                                  name=f"torchgwas-sigshard-{index}", daemon=True)
        thread.start()
        workers.append(thread)
    remaining = len(workers)
    try:
        while remaining:
            if failures:
                raise failures[0]
            try:
                item = results.get(timeout=0.05)
            except queue.Empty:
                continue
            if item is done:
                remaining -= 1
            else:
                yield item
    finally:
        stop.set()
        for thread in workers:
            thread.join()
    if failures:
        raise failures[0]

