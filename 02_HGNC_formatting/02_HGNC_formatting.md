02_miRNA_analysis: HGNC formatting
================
D. Ryazantsev, E. Sharova
2026-04-21

- [02 Анализ миРНК](#02-анализ-мирнк)
- [Объединение всех генов миРНК из чипа с соответствующими им генами из
  аннотации.](#объединение-всех-генов-мирнк-из-чипа-с-соответствующими-им-генами-из-аннотации)
  - [Какие есть различия в названиях:](#какие-есть-различия-в-названиях)
    - [Объединяем в общую таблицу](#объединяем-в-общую-таблицу)
- [Итоговая таблица](#итоговая-таблица)
- [Пересечем mimat gff3](#пересечем-mimat-gff3)
  - [Дэги из статьи](#дэги-из-статьи)
- [Делаем бедки из полученной таблицы
  миРНК](#делаем-бедки-из-полученной-таблицы-мирнк)
- [Поиск таргетов miRNA](#поиск-таргетов-mirna)

# 02 Анализ миРНК

Импортируем предыдущие данные и gtf

**Части анализа:**

1.  Получить пересечения всех генов из чипа с генами, по стренду и
    против
    1.  В границах гена: по стренду и против

    2.  На небольшом расстоянии от гена
2.  Посмотреть корреляцию с транскриптомами
3.  Смотреть таргеты
4.  Пересечь все с дифэкспрессированными

Импорт полученных ранее данных

``` r
# импортируем список дифэкспрессированных миРНК
diffexpressed_mirna_HGNC <- readRDS("/data7a/bio/human_genomics/fuchs_dystrophy/nanostring/analisys/RDS/diffexpressed_mirna_HGNC.RDS")
chip_names <- readRDS("/data7a/bio/human_genomics/fuchs_dystrophy/nanostring/analisys/RDS/chip_names.RDS")
norm_data <- readRDS("/data7a/bio/human_genomics/fuchs_dystrophy/nanostring/analisys/RDS/norm_data.RDS")
#annot <- as.data.frame(rtracklayer::import("/data7a/bio/human_genomics/shared/maps/GDC.h38.d1.vd1/gencode.v43.chr_patch_hapl_scaff.annotation.dic.gtf","gtf"))})

#saveRDS(annot, file="/data7a/bio/human_genomics/fuchs_dystrophy/nanostring/analisys/RDS/gencode.v43.chr_patch_hapl_scaff.annotation.dic.gtf.RDS")
# 37 vs 16 секунд, в 2.5 раза
annot <- readRDS("/data7a/bio/human_genomics/fuchs_dystrophy/nanostring/analisys/RDS/gencode.v43.chr_patch_hapl_scaff.annotation.dic.gtf.RDS")
#DT::datatable(chip_names, caption = "Таблица с микроРНК в микрочипе")
```

# Объединение всех генов миРНК из чипа с соответствующими им генами из аннотации.

Сначала делаем универсальную таблицу из аннотации, с пересекаемыми
именами для df со всеми МИР в чипе, потом пересекаем их, получая
финальную таблицу миРНК в чипе с колонками с другими ID из аннотации.

## Какие есть различия в названиях:

- y - решено

- res - ресёрч, принятие решения

- dev - решение в разработке.  
  mirbase mature ID (MIMAT) и pre ID (MI).

<table style="width:99%;">
<colgroup>
<col style="width: 6%" />
<col style="width: 20%" />
<col style="width: 42%" />
<col style="width: 30%" />
</colgroup>
<thead>
<tr class="header">
<th>Fix</th>
<th>chip_names (A)</th>
<th>annot (HGNC) (B)</th>
<th>Comment/Решение</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td>y</td>
<td>LET7A</td>
<td>MIRLET7A</td>
<td>исправить в <strong>(А)</strong></td>
</tr>
<tr class="even">
<td>y(29.10)</td>
<td>MIRLET7A<br />
MIRLET7F</td>
<td>MIRLET7A1,2,3</td>
<td>123 это шпильки, они идут по одному MIMAT добавить до записи, с
MIRLET7F тоже<strong><br />
</strong>Добавил вместе с a123 <strong>(B)</strong></td>
</tr>
<tr class="odd">
<td>res?</td>
<td><a href="https://mirbase.org/mature/MIMAT0004803">MIR548A
5P</a><br />
<a href="https://mirbase.org/mature/MIMAT0003251">MIR548A 3P</a></td>
<td>MIR548A1<br />
MIR548A2</td>
<td>в 3P шпильки (123), 5P только 3 (та же, что и в 5P). Пока не знаю,
чд</td>
</tr>
<tr class="even">
<td>y</td>
<td>MIR101</td>
<td>MIR101-[1,2,3]</td>
<td>суффиксы вынесены в отд колонку в <strong>(B)</strong>
<sup>1A*</sup></td>
</tr>
<tr class="odd">
<td><p>y(22.10)</p>
<p>(отмена 29.10)</p></td>
<td>MIR128-1</td>
<td>MIR128-1</td>
<td>в А сразу с -1, а в В мы вынесли -1 в отд колонку.<br />
Наверное нужно -1 в <strong>(А)</strong> вынести тоже в суфф.<br />
<u>UPD: буду метчить эти колонки во 2й подход.</u></td>
</tr>
<tr class="even">
<td>y(29.10)</td>
<td><ul>
<li><p>MIR103A<br />
то же</p></li>
<li><p>MIR19B</p></li>
</ul></td>
<td><ul>
<li><p>MIR103A<br />
A1 A2 A3</p></li>
<li><p>MIR19B1, MIR19B2</p></li>
</ul></td>
<td>Добавить все шпильки (куда?) (они без “-”).<br />
B1 и B2 с разных хромосом.<br />
Реш: вынести <strong>(В)</strong> в отд колонку энивей,
регуляркой.<br />
Лучше: вынести их в отд таблицу, только в ней всё поправить, и потом
rbind.</td>
</tr>
<tr class="odd">
<td>y(31.10)</td>
<td>MIR92A-1</td>
<td>MIR92A1</td>
<td>Испр в <strong>(А)</strong>, смержить отдельно по ним, потом
добавить.</td>
</tr>
<tr class="even">
<td>y(18.10)</td>
<td>MIR11852<br />
(дб 1185-2)</td>
<td>NA<br />
(непр имя)</td>
<td>исправить имя в <strong>(А)</strong> на 85-2 (2 знака -)</td>
</tr>
<tr class="odd">
<td>y(18.10)</td>
<td>MIR106A+MIR17</td>
<td>MIR106A и MIR17 or MIR106A</td>
<td>mimat у них один и выводит 1ю миР.<br />
Решение: вывести 1 в base_name,<br />
остальные <del>в суффикс</del> убрать.<br />
Возможно не оч решение, но пока ок</td>
</tr>
<tr class="even">
<td>y(29.10)</td>
<td>MIRLET7A</td>
<td>MIRLET7A[123]</td>
<td>MI0000062 - соотв 7A3</td>
</tr>
<tr class="odd">
<td>y(29.10)</td>
<td><a href="https://mirbase.org/mature/MIMAT0005905">MIR1254</a><br />
MIR4461<br />
MIR4532<br />
MIR4792<br />
MIR566</td>
<td>NA<br />
<a
href="https://www.genenames.org/data/gene-symbol-report/#!/hgnc_id/HGNC:41770">This
record has been withdrawn by HGNC</a></td>
<td>нет MIR имени, по mimat ищется.<br />
Шпильки -1 -2 отозваны из mirbase, HGNC.<br />
Решение - удалить из <strong>(А)</strong> /не добавлять.</td>
</tr>
</tbody>
</table>

> - 1A 🏃  
>   Решить, названия с суффиксами -1 -2, как объединяются с названием,
>   писать 1 в каждый -n или выбрать только 1 -n?
>
> - насчёт стрендов? - отметил их в (A).
>
> - нельзя добавить в таблицу name\*, тк там могут быть куча левых букв,
>   нужно прописать цифры в re \[1-9\].
>
> - Мнение насчёт шпилек (-1 -2 и 1,2): зонды в чипе идут на просто
>   *mature ID (MIMAT)* последовательность миРНК и мы не можем
>   специфично сказать, с какого места в геноме эта миРНК
>   экспрессируется, поэтому единственное, что мы можем сделать -
>   добавить все эти шпильки *pre ID (MI)* в общую таблицу.  
>   Эти шпильки практически идентичны по длине (±5bp), но могут
>   экспрессироваться с разных хромосом. Но мы ничего тут поделать не
>   можем.  
>   **Варианты добавляния их:**
>
>   - Добавлять из (А) их как регулярки вида “MIR103A.”, чтобы добавлять
>     посл символ, если это цифра или любой 1 символ или “-символ”.
>
>   - Поделить в (В) ещё раз суффиксы без “-”, чтобы оставить base_name
>     и suffix. И тогда по base_name оставить вхождения в (В).

``` r
# сделаем подтаблицу из аннотации, где будут только мирнк гены
# возможно  мы  упустим так пару генов,которые не МИР но проверим это позже
# разобьём сложные имена из аннотации по столбцам, чтобы объединить их
# далее работаем с этим df
annot_mirdf <- annot[annot$gene_type == 'miRNA', ]



# (B) MIR123-2 --> MIR123, 2
# Создаем новые колонки с разделением на base и suffix
# Разделяем gene_name на две части и сохраняем в новые колонки
annot_mirdf <- separate(annot_mirdf, gene_name, into = c("name_base", "name_suffix"), sep = "-", remove = F)


# (B) MIR19B1 --> MIR19B, 1
# (B) MIRLET7A1 --> MIRLET7A, 1 (заодно)
# (B) MIR548AA1 --> MIR548AA, 1 (заодно)
# выделим в отдельную подтаблицу, в ней исправим регуляркой 
# потом соединим их обратно 
# выделяем в отд таблицу
filtered_annot_mirdf <- annot_mirdf %>%
  filter(grepl("MIR(LET)*[0-9]*[A-Z]{1,2}[0-9]", gene_name))

# удаление выбранных строк из основной таблицы
annot_mirdf <- annot_mirdf %>%
  filter(!grepl("MIR(LET)*[0-9]*[A-Z]{1,2}[0-9]", gene_name))

# Разделяем по колонкам
filtered_annot_mirdf <- filtered_annot_mirdf %>%
  extract(gene_name, into = c("name_base", "name_suffix"), regex = "([A-Z0-9]+[A-Z]{1,2})([0-9])", remove = F)

# соединяем таблицы обратно
annot_mirdf <- bind_rows(annot_mirdf, filtered_annot_mirdf)
```

Сначала посмотрим какие гены не пересекаются с аннотацией и почему

``` r
# проверяем насколько пересечётся аннотация с генами из чипа
chip_diff <- annot_mirdf[annot_mirdf$name_base %in% chip_names$HGNC_name, ] #1971 657
not_intersected_in_chip <- setdiff(chip_names$HGNC_name, chip_diff$name_base)
not_intersected_in_chip
```

    ##  [1] "MIR1185-1" "MIR1185-2" "MIR1254"   "MIR128-1"  "MIR128-2"  "MIR129-2" 
    ##  [7] "MIR181A2"  "MIR181B2"  "MIR219A1"  "MIR219A2"  "MIR376A2"  "MIR4461"  
    ## [13] "MIR450A1"  "MIR450A2"  "MIR4532"   "MIR4792"   "MIR509-3"  "MIR566"   
    ## [19] "MIR92A1"

``` r
length(not_intersected_in_chip)
```

    ## [1] 19

``` r
rm(chip_diff)
#18.10 - 72
#18.10 - 54
#22.10 - 40
#29.10 - 22 (a123)
#29.10 - 20 (mirlet a123)
#29.10 - 15 (удалил отозв-е из анота)
#29.10 - 14 ( )
# осталось 14 шт, смержим их отдельно по полному имени, а не по name_base
```

Далее объединяем подготовленную таблицу из аннотации с именами из чипов.
Нужно сметчить по имени “MIR123” ака HGNC_name. Ещё можно убрать совсем
пустые колонки (или не нужно?).

Варианты объединения:

- Выделять ли из аннотации лишь нужные вхождения

- подтянуть в аннотацию всю таблицу из чипов?  
  Сделал та, колонка name_base удалилась, теперь только hgnc. Пока
  оставил так, но надо спросить и пофиксить, сделав как в 1.

  с чертой мержнуть в отд таблицу, а потом 2 rbind

### Объединяем в общую таблицу

``` r
# мерджим по полному названию gene_name
merge_names <- chip_names %>% 
  filter(grepl("MIR[0-9]*-[0-9]|MIR[0-9]*[A-Z]{1,2}[0-9]", HGNC_name))

merge_df <- merge(merge_names, annot_mirdf, 
                            by.x = "HGNC_name", by.y = "gene_name", all.y = F, all.x = F)

# удаляем эти вхождения из аннотации, чтобы не вытянуть лишние след шагом
annot_mirdf <- anti_join(annot_mirdf, merge_df, by = c("gene_name" = "HGNC_name"))
# 5637 - 42 == 5595
```

``` r
# комбинировать в 2 этапа сначала с "-" потом без
# мержим обычные HGNC_name+name_base
annot_chipdf <- merge(chip_names, annot_mirdf, 
                            by.x = "HGNC_name", by.y = "name_base", all.y = F, all.x = F)
# 2709 

# гены через - подтянулись просто доп строчками
# по gene_name подтягивается по всем, что неправильно тк там все подрят
# annot_chipdf <- merge(chip_names, annot_mirdf, 
#                             by.x = "HGNC_name", by.y = "gene_name", all.y = F, all.x = F)

# соединяем 1 и 2 мерж
annot_chipdf <- bind_rows(annot_chipdf, merge_df)
# 2751 
```

Проверяем как получилось

    1
    MIR873 5р 3р -> все гуд отдельро под каждый стренд запись
    MIR1283 -> MIR1283 -2 отдельтно под стренды
    MIR103A -> MIR103A1,2 отдельно под каждый А123
    MIR129-2 -> подтянулись все по имени 129, то есть лишнее. у записей -2 другой стренд и другой mimat-id
    -подправил все и все ок

> Синхронизировать по названию выходит проблематично, поэтому можно было
> найти таблицу соответствия по MIMAT id. Нашёл файл формат embl, но не
> устанавливается пакет, читающий его. На будующее `gbRecord` — пакет
> для импорта EMBL формата, тул для переименования миРНК —
> `miRNAmeConverter`.

# Итоговая таблица

Мы подтянули данные из аннотации для каждого вхождения микроРНК (и
подтипов) в микрочипе, расширив изначал

Таблица, где только гены. Объединил колонки name_base, gene_name в одну

Подсчёт уникальных и повторяющихся вхождений.

# Пересечем mimat gff3

Скачал файл hsa.gff3 с [mirbase](https://www.mirbase.org/download/). Там
уже есть MIMAT ID-шники, сметчим по ним. Не помню, почему раньше не
находил его. Помню что скачал оттуда большую таблицу и искал там что-то
через libre office и не нашёл.

### Дэги из статьи

Пересечём дэги из статьи с аннотацией. Проблема: много генов не
пересекаются, тк они указаны как mature miRNA (hsa-mi<u>**R**</u>-26b),
но суффикса -5/3 у большинство них нет, поэтому много теряется. В итоге
получается 80/311 (25%), оставим пока как есть. UPD 11.08.25 после
получения полных названий пересеклись всё и даже больше ( \_1, 2, 3)

Вообще я пересекал bedtools только HGNC перевод, думаю стоит пересечь и
mimat табличку., в этот раз добавив ей оба mimat id

``` r
matthaei_degs <- readRDS("/data7a/bio/human_genomics/fuchs_dystrophy/nanostring/analisys/RDS/matthaei_degs.RDS")
matthei_gff <- left_join(matthaei_degs, mimat_table, by = c("miRNA" = "Name")) %>% 
  select(-c(source, score, phase)) %>% filter(!is.na(seqnames)) %>% 
  distinct() %>% mutate(uniq_pos=paste(seqnames, start, end, sep="_" ))

saveRDS(matthei_gff, file="/data7a/bio/human_genomics/fuchs_dystrophy/nanostring/analisys/RDS/matthei_gff.RDS")
```

``` r
chip_mimat <-  merge(chip_names, mimat_table, by.y="Alias", by.x="MIMAT_ID") %>% 
  select(-c(source, score, phase)) %>% 
  distinct() %>% mutate(uniq_pos=paste(seqnames, start, end, sep="_" )) %>% 
  
  group_by(MIMAT_ID) %>% 
  mutate(dups_mimat = n()) %>% 
  ungroup() %>% 
  relocate(dups_mimat, .before = everything()) %>% 
  
  group_by(miRNA) %>% 
  mutate(dups_miRNA = n()) %>% 
  ungroup() %>% 
  relocate(dups_miRNA, .before = everything())

# не пересеклись
chip_names[!(chip_names$MIMAT_ID %in% mimat_table$Alias), ]
```

    ##            miRNA HGNC_name mir_strand     MIMAT_ID
    ## 55  hsa-miR-1254   MIR1254         5P MIMAT0005905
    ## 207 hsa-miR-1973   MIR1973         5P MIMAT0009448
    ## 408 hsa-miR-378g   MIR378G         5P MIMAT0018937
    ## 451 hsa-miR-4461   MIR4461         5P MIMAT0018983
    ## 469 hsa-miR-4532   MIR4532         5P MIMAT0019071
    ## 482 hsa-miR-4792   MIR4792         5P MIMAT0019964
    ## 639  hsa-miR-566    MIR566         5P MIMAT0003230

``` r
# одинаковы ли повторения?
unique(chip_mimat$dups_miRNA == chip_mimat$dups_miRNA)
```

    ## [1] TRUE

``` r
#(chip_mimat[(chip_mimat$Name %in% targets_table$miRNA), ])
# одинаковы ли имена
#DT::datatable(chip_mimat[!chip_mimat$miRNA == chip_mimat$Name,] %>% select(miRNA, Name), caption = 'Сравнение имен из чипа и MIMAT gff3 в получишейся таблице')
```

Имена практически совпадают, из не совпадающих - мир с плюсами в
названии, иногда есть приписка -5\3P. Итог - мержить надо по названиям
из столбца Name(Mimat gff3).

Получившаяся таблица

Не пересеклись с mimat таблицой 2 шт hsa-miR-1973, hsa-miR-378g.

тоже сделаем статистику по повторяющимся вхождениям. По колонкам Name и
miRNA одинаково. Повторения мимат айдишников и названий абсолютно
идентичны.

``` r
uniq_chip_mimat <- full_join(
  chip_mimat %>% dplyr::count(dups_miRNA) %>% rename(count = dups_miRNA, dups_miRNA_n = n),
  chip_mimat %>% dplyr::count(dups_mimat) %>% rename(count = dups_mimat, dups_mimat_n = n),
  by = "count"
) %>% replace_na(list(dups_ENSG_n = 0, dups_mimat_n = 0))


# почему-то в for article rmd выдавал ошибку на knit-е
# хотя отдельно каждый из скриптов сорсился
saveRDS(uniq_chip_mimat, file="/data7a/bio/human_genomics/fuchs_dystrophy/nanostring/analisys/RDS/uniq_chip_mimat.RDS")


#DT::datatable(uniq_chip_mimat, caption = "Сколько раз повторяются запись МИР и MIMAT, пересеч с gff списком MIMAT")
```

# Делаем бедки из полученной таблицы миРНК

Из полученной таблицы `annot_chipdf` выделим только гены, получим из них
bed-файл. Нуж перевести в bed 0-based формат(в GFF - 1-based), поэтому
start - 1

``` r
# make bed

annot_chipdf_bed <- annot_chipdf_genes[, c("seqnames", "start", "end", "gene_id","name_base","strand","mir_strand")]
chip_mimat_bed <- chip_mimat[, c("seqnames", "start", "end", "miRNA", "HGNC_name","strand", "MIMAT_ID", "ID","Derives_from", "mir_strand")] %>% mutate(start = start - 1)

matthei_mimat_bed <- matthei_gff[, c("seqnames", "start", "end", "miRNA", "type", "strand.y", "Alias", "ID", "Derives_from",   "uniq_pos")] %>% mutate(Alias = as.character(Alias)) %>% mutate(start = start - 1)

write.table(annot_chipdf_bed, file="/data7a/bio/human_genomics/fuchs_dystrophy/nanostring/analisys/download/annot_chipdf.bed", quote = F, sep="\t", col.names =F, row.names = F )
write.table(chip_mimat_bed, file="/data7a/bio/human_genomics/fuchs_dystrophy/nanostring/analisys/download/chip_mimat.bed", quote = F, sep="\t", col.names =F, row.names = F )

write.table(matthei_mimat_bed, file="/data7a/bio/human_genomics/fuchs_dystrophy/nanostring/analisys/download/matthaei_mimat.bed", quote = F, sep="\t", col.names =F, row.names = F )

#DT::datatable(annot_chipdf_bed, caption = "bed файл")
```

``` r
annot_chipdf_bed <- annot_chipdf_genes[, c("seqnames", "start", "end", "gene_id","name_base","strand","mir_strand")]


write.table(annot_chipdf_bed, file="/data7a/bio/human_genomics/fuchs_dystrophy/nanostring/analisys/download/annot_chipdf_bed.bed", quote = F, sep="\t", col.names =F, row.names = F )
```

Сравнить стренд 5Р/3Р и +/- в аннотации – они различаются. Планировалось
сделать если 5p, то стренд как в аннотации, если 3p, то реверснуть как в
аннотации. (пока не сделано)

Получаем аннотационный файл, из которого удалим чиповые миРНК. Удалить
абсолютно все быстро и легко не вышло, пока оставим так. Это можно будет
подправить на постобработке, если понадобится.

Мирнк внутри генов будет коэкспрессироваться вместе с ними и
положительно коррелировать с утровнем экспрессии, которую мы посмотрим в
соответствующих транскриптомах.

``` r
# отбираем из аннотационного ф-ла гены миРНК и аннотацию генов без миРНК
# все миРНК
# an_mirna <- annot[annot$gene_type=="miRNA" & annot$type=="gene", ]
# an_mirna <- an_mirna[, c("seqnames", "start","end", "gene_id","gene_name","strand","gene_type")]
# write.table(an_mirna, file="/data7a/bio/human_genomics/fuchs_dystrophy/nanostring/analisys/download/an_mirna.bed", quote = F, sep="\t", col.names =F, row.names = F )

an_wo_mirna <- anti_join(annot, annot_chipdf, by = c("gene_name" = "HGNC_name"))

an_wo_mirna <- an_wo_mirna[, c("seqnames", "start", "end", "gene_id","gene_name", "strand","type", "gene_type")]

write.table(an_wo_mirna, file="/data7a/bio/human_genomics/fuchs_dystrophy/nanostring/analisys/download/an_wo_mirna.bed", quote = F, sep="\t", col.names =F, row.names = F)
```

**Получаем bed пересечения:**

- Внутри генов - получаем миртроны, коэкспрессированные с генов

  - по и против стренда

- Около генов на расстоянии 10.000bp  
  Увидим открытые участки хроматина, которые

  - по и против стренда

``` bash
#!/bin/bash
cd "/data7a/bio/human_genomics/fuchs_dystrophy/nanostring/analisys/download/"
eval "$(conda shell.bash hook)"
conda activate bedtools

# без учета strand
bedtools intersect -wb -a an_wo_mirna.bed -b annot_chipdf_genes.bed > annot_chipdf_genes_s_no.bed
# -s strand-specific (оба на + или оба на -)
bedtools intersect -wb -s -a an_wo_mirna.bed -b annot_chipdf_genes.bed > annot_chipdf_genes_s_y.bed
# -S strand-opposite (один на +, другой на -)
bedtools intersect -wb -S -a an_wo_mirna.bed -b annot_chipdf_genes.bed > annot_chipdf_genes_s_reverse.bed

# около генов
# без учета strand
bedtools window -u -a an_wo_mirna.bed -b annot_chipdf_genes.bed -w 10000 > neargene_annot_chipdf_genes_s_no.bed 
# -sm same strand strand-specific (оба на + или оба на -)
bedtools window -u -sm -a an_wo_mirna.bed -b annot_chipdf_genes.bed -w 10000 > neargene_annot_chipdf_genes_s_y.bed 

# новое
bedtools intersect -wb -s -a an_wo_mirna.bed -b chip_mimat.bed > chip_mimat_s_y.bed
bedtools intersect -wb -s -a an_wo_mirna.bed -b matthei_mimat.bed > matthei_mimat_s_y.bed
```

Забираем получившиеся пересечения для дальнейшего анализа.

``` r
genes_s_y <- read_delim("/data7a/bio/human_genomics/fuchs_dystrophy/nanostring/analisys/download/annot_chipdf_genes_s_y.bed", delim="\t", col_names=F)
genes_s_n <- read_delim("/data7a/bio/human_genomics/fuchs_dystrophy/nanostring/analisys/download/annot_chipdf_genes_s_no.bed", delim="\t", col_names=F)
neargene_s_y <- read_delim("/data7a/bio/human_genomics/fuchs_dystrophy/nanostring/analisys/download/neargene_annot_chipdf_genes_s_y.bed", delim="\t", col_names=F)
neargene_s_n <- read_delim("/data7a/bio/human_genomics/fuchs_dystrophy/nanostring/analisys/download/neargene_annot_chipdf_genes_s_no.bed", delim="\t", col_names=F)
```

Отфильтруем миРНК из аннотации, которые мапятся сами на себя по названию
хромосоме и координатам, оставляем уникальные вхождения. Вхождения вида
MIR-HG (host genes) не отфильтровывались, но их дб в обоих таблицах
одинаково, так что пока пусть будут.

``` r
genes_s_y <- genes_s_y %>%
  filter(!(X1 == X9 & X2 == X10 & X3 == X11), X7 == 'gene') %>% distinct()
genes_s_n <- genes_s_n %>%
  filter(!(X1 == X9 & X2 == X10 & X3 == X11), X7 == 'gene') %>% distinct()
```

Делаем итоговую таблицу, сколько генов мы получили

``` r
table_genes <- tibble(
"По стрэнду" =  length(unique((genes_s_y$X5))),
"Против стрэнда" =  length(unique((genes_s_n$X5)))) 
table_genes
```

    ## # A tibble: 1 × 2
    ##   `По стрэнду` `Против стрэнда`
    ##          <int>            <int>
    ## 1          432              566

Пока берём паузу на этом этапе.

# Поиск таргетов miRNA

<u>**UPD**</u>**:** для дальнейшего анализа этот этап опустили. Брали
таблицу с таргетами целиком, без фильтраций по значимости.

Импортируем файл с таргетами из miRTarBase Release 8.0 ([нашли его в
нашей
статье](https://pmc.ncbi.nlm.nih.gov/articles/PMC10056610/#B54-life-13-00659)).
[Ссылка на
сайт](https://mirtarbase.cuhk.edu.cn/~miRTarBase/miRTarBase_2019/php/download.php?ver=8.0&opt=show)

Есть новая база версии 10
[ссылка](https://awi.cuhk.edu.cn/~miRTarBase/miRTarBase_2025/php/download.php)

``` r
targets_table <- read_excel("/data7a/bio/human_genomics/fuchs_dystrophy/nanostring/analisys/upload/hsa_MTI_target_table.xlsx")
```

``` r
# фильтруем таблицу с таргетами, получаем более значимые
targets_table <- targets_table %>% filter(`Support Type` == "Functional MTI" | `Support Type` == "Functional MTI (Weak)" )

# имена из чипа пересекаются на 777/798 с таргетами
length(intersect(targets_table$miRNA, chip_names$miRNA))
not_intersected_in_target <- setdiff(chip_names$miRNA, targets_table$miRNA)
# значения вида mir+mir (только 2шт) не пересекаются, поэтому разделим их на 2 строчки
chip_names <- chip_names %>%
  separate_rows(miRNA, sep = "\\+")

# готово, все пересекаются not_intersected_in_target = 0
not_intersected_in_target <- setdiff(chip_names$miRNA, targets_table$miRNA)



# получаем таблицу с таргетами по миРНК чипам
chip_targets <- merge(chip_names, targets_table, 
                            by.x = "miRNA", by.y = "miRNA", all.y = F, all.x = F)
```

Ищем таргеты для дифэкспрессированных миРНК.

``` r
targets_diffexpr <- merge(diffexpressed_mirna_HGNC, targets_table, 
                            by.x = "miRNA", by.y = "miRNA", all.y = F, all.x = F)
```

❗ В списке с u-тестом есть мирнк, подходящая по таргетам на TCF4, в
лимме и т-тесту пересечений нет. <u>Ответ</u>: Лучше так не делать, TCF4
оч длинный ген и нет смысла его фокусить по миРНКе.

    "hsa-miR-181a-5p" "hsa-miR-519d-3p"

Ищем координаты миРНК в генах

Сначала считаем на всех миРНК, потом на дифэкспрессированных.

Из аннотационного файла получить bed-ки, в 1 вытащить все миРНК с их
координатами, в 2 мы сначала убираем все миРНК, потом оставляем
координаты со strand и id гена.

Всех микрнк чипа, убрать их из аннотации

Про дифэкспрессирующиеся мирнк: не все вытаскиваются. Их нет в
аннотации: MIR128 “MIR1283

1.  MIR128 просто нет, есть MIR128-1, MIR128-2. В таблице идет как
    hsa-miR-128-3p.  
    \* 128-1, 128-2 значит, что один сиквенс с разных локусов.  
    Итог: Лучше добавить обе записи (128-1, 128-2), тк обе
    [ссылаются](https://www.genecards.org/cgi-bin/carddisp.pl?gene=MIR128-1&keywords=MIR128-1)
    на название hsa-miR-128-3p. Это
    [шпильки](https://mirbase.org/mature/MIMAT0000424) 128

2.  1283-1 -2 это шпильки 1283, добавить обе

    Добавляем записи 128-1/2 вручную: MIR128-1, MIR128-2, MIR1283-1,
    MIR1283-2

**Дальнейший анализ в файле 03\_**
