import os
import copy
import pandas as pd

from gcsnap.configuration import Configuration
from gcsnap.rich_console import RichConsole
from gcsnap.genomic_context import GenomicContext
from gcsnap.annotations.operons import Operons

class DeNovo:
    """
    Assign protein families to flanking genes from an external annotation file,
    bypassing sequence clustering.

    Each unique annotation_id in the file defines one family (numbered 1, 2, 3, ...,
    ordered by decreasing frequency in the annotation file). Proteins absent from
    the annotation file are placed in family 0 (unannotated). Pseudogenes retain
    family -1. The target gene always receives max_family + 1.

    Protein names in the copied syntenies are replaced by annotation_name wherever a
    match exists (pseudogenes and the target gene are never renamed).

    All output files are written to out_label (a fresh, independent directory).

    Attributes:
        config (Configuration): Configuration object.
        annotation_file (str): Path to the TSV annotation file.
        out_label (str): Output directory for the annotated genomic context.
        annotated_gc (GenomicContext): The new annotated GenomicContext.
        cds_to_annotation (dict): Mapping cds_code -> (annotation_id, annotation_name).
        annotation_id_to_family (dict): Mapping annotation_id -> int family number.
        console (RichConsole): Console for printing messages.
    """

    def __init__(
        self,
        config: Configuration,
        gc: GenomicContext,
        annotation_file: str,
        out_label: str,
    ) -> None:
        """
        Initialize the DeNovo object.

        Args:
            config (Configuration): Configuration object.
            gc (GenomicContext): Source GenomicContext (not modified).
            annotation_file (str): Path to a TSV file with columns:
                cds_code, annotation_id, annotation_name, confidence.
            out_label (str): Directory where annotated output files will be written.
        """
        self.config = config
        self.out_label = out_label
        self.annotation_file = annotation_file
        self.console = RichConsole()

        os.makedirs(out_label, exist_ok=True)

        # Create a fresh GenomicContext pointing to out_label, then populate
        # it with a deep copy of the source syntenies and taxonomy so that the
        # source gc is never modified and all file paths point to out_label.
        self.annotated_gc = GenomicContext(config, out_label)
        self.annotated_gc.syntenies = copy.deepcopy(gc.syntenies)
        self.annotated_gc.taxonomy = copy.deepcopy(gc.taxonomy)
        self.annotated_gc.curr_targets = list(gc.curr_targets)

        # Read annotation file
        annotations = pd.read_csv(annotation_file, sep='\t')

        # cds_code -> (annotation_id, annotation_name)
        self.cds_to_annotation = {
            row['cds_code']: (row['annotation_id'], row['annotation_name'])
            for _, row in annotations.iterrows()
        }

        # annotation_id -> family number (1-based, most frequent first)
        id_counts = annotations['annotation_id'].value_counts()
        self.annotation_id_to_family = {
            aid: i + 1 for i, aid in enumerate(id_counts.index)
        }

    # ------------------------------------------------------------------
    # Public interface
    # ------------------------------------------------------------------

    def run(self) -> None:
        """
        Full pipeline:
            1. Assign families and update protein names in syntenies.
            2. Build and write the families summary.
            3. Re-run operon clustering on the annotated genomic context.
            4. Write operon summary and final syntenies JSON.
        """
        with self.console.status('Assigning de novo families from annotation file'):
            self._assign_families()

        # Build gc.families from the updated syntenies (reset first to avoid duplicates)
        self.annotated_gc.families = {}
        self.annotated_gc.create_families_summary()
        self.annotated_gc.write_families_summary_to_txt()
        self.annotated_gc.write_families_to_json()

    def get_annotated_gc(self) -> GenomicContext:
        """Return the fully annotated GenomicContext."""
        return self.annotated_gc

    # ------------------------------------------------------------------
    # Internal helpers
    # ------------------------------------------------------------------

    def _assign_families(self) -> None:
        """
        Walk every flanking gene in the copied syntenies and:
        - Rename the protein to annotation_name if present in the annotation file.
        - Assign the family number that corresponds to its annotation_id.
        - Leave pseudogenes at family -1 (names untouched).
        - Give the target gene family max_family + 1 (name untouched).
        - Assign family 0 to all remaining (unannotated) proteins.

        The priority mirrors Families.assign_families():
            pseudogene check runs first, target-gene check overrides it.
        """
        max_family = (
            max(self.annotation_id_to_family.values())
            if self.annotation_id_to_family
            else 0
        )

        for target_key, target_data in self.annotated_gc.syntenies.items():
            flanking = target_data['flanking_genes']
            families_list = []

            # The actual CDS code of the target protein in the assembly
            # (may differ from target_key for some providers).
            target_cds = target_data.get('assembly_metadata', {}).get(
                'target_cds', target_key
            )

            for i, cds_code in enumerate(flanking['cds_codes']):
                name = flanking['names'][i]

                # Default: unannotated
                family = 0

                # 1. Pseudogene
                if name == 'pseudogene':
                    family = -1

                # 2. Known annotation (skipped for pseudogenes)
                if cds_code in self.cds_to_annotation and name != 'pseudogene':
                    ann_id, ann_name = self.cds_to_annotation[cds_code]
                    family = self.annotation_id_to_family[ann_id]
                    flanking['names'][i] = ann_name
                elif name != 'pseudogene':
                    # Not in annotation file and not a pseudogene: rename to unannotated
                    flanking['names'][i] = 'Unannotated protein'

                # 3. Target gene always overrides (mirrors Families behaviour)
                if cds_code == target_key:
                    family = max_family + 1

                families_list.append(family)

                # Record target_family based on assembly_metadata.target_cds
                if cds_code == target_cds:
                    target_data['target_family'] = family

            flanking['families'] = families_list