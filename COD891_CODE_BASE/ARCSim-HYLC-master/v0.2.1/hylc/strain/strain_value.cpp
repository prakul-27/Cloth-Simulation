#include "strain.hpp"
#ifndef hylc_strain_II

using namespace hylc;
using namespace hylc::mathematica;

Vec6 hylc::mathematica::strain(const Vec18 &xloc, const Mat2x2 &invDm,
                               const Vec3 &niexists) {
  // define output
  Vec6 out;

  Real copt42  = invDm(0, 0);
  Real copt45  = xloc(0);
  Real copt47  = -copt45;
  Real copt49  = xloc(3);
  Real copt50  = copt47 + copt49;
  Real copt54  = copt42 * copt50;
  Real copt57  = invDm(1, 0);
  Real copt62  = xloc(6);
  Real copt64  = copt47 + copt62;
  Real copt65  = copt57 * copt64;
  Real copt68  = copt54 + copt65;
  Real copt70  = Power(copt68, 2);
  Real copt72  = xloc(1);
  Real copt73  = -copt72;
  Real copt77  = xloc(4);
  Real copt78  = copt73 + copt77;
  Real copt81  = copt42 * copt78;
  Real copt86  = xloc(7);
  Real copt189 = copt73 + copt86;
  Real copt195 = copt189 * copt57;
  Real copt207 = copt195 + copt81;
  Real copt236 = Power(copt207, 2);
  Real copt241 = xloc(2);
  Real copt244 = -copt241;
  Real copt251 = xloc(5);
  Real copt270 = copt244 + copt251;
  Real copt297 = copt270 * copt42;
  Real copt298 = xloc(8);
  Real copt299 = copt244 + copt298;
  Real copt306 = copt299 * copt57;
  Real copt311 = copt297 + copt306;
  Real copt332 = Power(copt311, 2);
  Real copt333 = copt236 + copt332 + copt70;
  Real copt334 = Sqrt(copt333);
  Real copt335 = 1 / copt334;
  Real copt336 = invDm(0, 1);
  Real copt355 = invDm(1, 1);
  Real copt350 = copt336 * copt50;
  Real copt356 = copt355 * copt64;
  Real copt357 = copt350 + copt356;
  Real copt393 = Power(copt357, 2);
  Real copt361 = copt336 * copt78;
  Real copt365 = copt189 * copt355;
  Real copt366 = copt361 + copt365;
  Real copt394 = Power(copt366, 2);
  Real copt368 = copt270 * copt336;
  Real copt382 = copt299 * copt355;
  Real copt387 = copt368 + copt382;
  Real copt395 = Power(copt387, 2);
  Real copt396 = copt393 + copt394 + copt395;
  Real copt397 = Sqrt(copt396);
  Real copt398 = 1 / copt397;
  Real copt403 = -copt49;
  Real copt404 = copt403 + copt45;
  Real copt405 = Power(copt404, 2);
  Real copt406 = -copt77;
  Real copt407 = copt406 + copt72;
  Real copt408 = Power(copt407, 2);
  Real copt409 = -copt251;
  Real copt410 = copt241 + copt409;
  Real copt412 = Power(copt410, 2);
  Real copt414 = copt405 + copt408 + copt412;
  Real copt415 = Sqrt(copt414);
  Real copt416 = 1 / copt415;
  Real copt417 = -copt62;
  Real copt418 = copt417 + copt45;
  Real copt419 = Power(copt418, 2);
  Real copt420 = -copt86;
  Real copt421 = copt420 + copt72;
  Real copt422 = Power(copt421, 2);
  Real copt423 = -copt298;
  Real copt424 = copt241 + copt423;
  Real copt426 = Power(copt424, 2);
  Real copt427 = copt419 + copt422 + copt426;
  Real copt428 = Sqrt(copt427);
  Real copt429 = 1 / copt428;
  Real copt430 = copt417 + copt49;
  Real copt431 = Power(copt430, 2);
  Real copt432 = copt420 + copt77;
  Real copt433 = Power(copt432, 2);
  Real copt436 = copt251 + copt423;
  Real copt437 = Power(copt436, 2);
  Real copt438 = copt431 + copt433 + copt437;
  Real copt439 = Sqrt(copt438);
  Real copt440 = 1 / copt439;
  Real copt441 = Power(copt42, 2);
  Real copt442 = Power(copt45, 2);
  Real copt443 = Power(copt72, 2);
  Real copt444 = Power(copt241, 2);
  Real copt446 = -2 * copt45 * copt49;
  Real copt447 = Power(copt49, 2);
  Real copt448 = -2 * copt72 * copt77;
  Real copt449 = Power(copt77, 2);
  Real copt450 = -2 * copt241 * copt251;
  Real copt451 = Power(copt251, 2);
  Real copt452 = copt442 + copt443 + copt444 + copt446 + copt447 + copt448 +
                 copt449 + copt450 + copt451;
  Real copt453 = copt441 * copt452;
  Real copt454 = -(copt241 * copt251);
  Real copt455 = copt49 * copt62;
  Real copt456 = copt49 + copt62;
  Real copt457 = -(copt45 * copt456);
  Real copt459 = copt77 * copt86;
  Real copt460 = copt77 + copt86;
  Real copt461 = -(copt460 * copt72);
  Real copt462 = -(copt241 * copt298);
  Real copt463 = copt251 * copt298;
  Real copt464 = copt442 + copt443 + copt444 + copt454 + copt455 + copt457 +
                 copt459 + copt461 + copt462 + copt463;
  Real copt466 = 2 * copt42 * copt464 * copt57;
  Real copt467 = Power(copt57, 2);
  Real copt469 = -2 * copt45 * copt62;
  Real copt470 = Power(copt62, 2);
  Real copt471 = -2 * copt72 * copt86;
  Real copt472 = Power(copt86, 2);
  Real copt473 = -2 * copt241 * copt298;
  Real copt474 = Power(copt298, 2);
  Real copt475 = copt442 + copt443 + copt444 + copt469 + copt470 + copt471 +
                 copt472 + copt473 + copt474;
  Real copt476 = copt467 * copt475;
  Real copt477 = copt453 + copt466 + copt476;
  Real copt478 = 1 / copt477;
  Real copt479 = copt430 * copt72;
  Real copt480 = copt62 * copt77;
  Real copt481 = -(copt49 * copt86);
  Real copt482 = copt406 + copt86;
  Real copt483 = copt45 * copt482;
  Real copt484 = copt479 + copt480 + copt481 + copt483;
  Real copt485 = Power(copt484, 2);
  Real copt486 = copt241 * copt430;
  Real copt487 = copt251 * copt62;
  Real copt488 = -(copt298 * copt49);
  Real copt489 = copt298 + copt409;
  Real copt491 = copt45 * copt489;
  Real copt492 = copt486 + copt487 + copt488 + copt491;
  Real copt493 = Power(copt492, 2);
  Real copt494 = copt241 * copt432;
  Real copt495 = copt251 * copt86;
  Real copt496 = -(copt298 * copt77);
  Real copt497 = copt489 * copt72;
  Real copt498 = copt494 + copt495 + copt496 + copt497;
  Real copt499 = Power(copt498, 2);
  Real copt500 = copt485 + copt493 + copt499;
  Real copt501 = Sqrt(copt500);
  Real copt502 = xloc(12);
  Real copt503 = -copt502;
  Real copt504 = copt49 + copt503;
  Real copt506 = copt417 + copt502;
  Real copt508 = copt403 + copt62;
  Real copt509 = xloc(13);
  Real copt515 = xloc(14);
  Real copt505 = copt504 * copt86;
  Real copt507 = copt506 * copt77;
  Real copt510 = copt508 * copt509;
  Real copt511 = copt505 + copt507 + copt510;
  Real copt557 = copt42 + copt57;
  Real copt558 = Power(copt557, 2);
  Real copt560 = Sqrt(copt452);
  Real copt561 = Sqrt(copt475);
  Real copt563 = -2 * copt49 * copt62;
  Real copt564 = -2 * copt77 * copt86;
  Real copt565 = -2 * copt251 * copt298;
  Real copt566 = copt447 + copt449 + copt451 + copt470 + copt472 + copt474 +
                 copt563 + copt564 + copt565;
  Real copt567 = Sqrt(copt566);
  Real copt568 = xloc(15);
  Real copt569 = -copt568;
  Real copt571 = copt569 + copt62;
  Real copt573 = copt47 + copt568;
  Real copt575 = xloc(16);
  Real copt581 = xloc(17);
  Real copt611 = -(copt62 * copt77);
  Real copt612 = copt432 * copt45;
  Real copt613 = copt49 * copt86;
  Real copt630 = xloc(9);
  Real copt631 = -copt630;
  Real copt637 = xloc(10);
  Real copt645 = xloc(11);
  Real copt528 = -(copt251 * copt45 * copt86);
  Real copt529 = copt298 * copt45 * copt77;
  Real copt632 = copt49 + copt631;
  Real copt544 = copt480 + copt481 + copt483;
  Real copt643 = copt403 + copt630;
  Real copt730 = copt500 * copt501;
  Real copt738 = copt333 * copt500;
  Real copt739 = Sqrt(copt738);
  Real copt740 = copt738 * copt739;
  Real copt741 = 1 / copt740;
  Real copt512 = copt484 * copt511;
  Real copt513 = copt298 * copt504;
  Real copt514 = copt251 * copt506;
  Real copt516 = copt508 * copt515;
  Real copt517 = copt513 + copt514 + copt516;
  Real copt518 = copt492 * copt517;
  Real copt519 = -copt509;
  Real copt520 = copt519 + copt77;
  Real copt521 = copt298 * copt520;
  Real copt522 = copt420 + copt509;
  Real copt523 = copt251 * copt522;
  Real copt524 = copt482 * copt515;
  Real copt525 = copt521 + copt523 + copt524;
  Real copt526 = copt498 * copt525;
  Real copt527 = copt512 + copt518 + copt526;
  Real copt530 = copt251 * copt502 * copt86;
  Real copt531 = -(copt298 * copt502 * copt77);
  Real copt532 = copt251 * copt45 * copt509;
  Real copt533 = -(copt251 * copt509 * copt62);
  Real copt534 = -(copt298 * copt45 * copt509);
  Real copt535 = copt298 * copt49 * copt509;
  Real copt536 = copt241 * copt511;
  Real copt545 = copt515 * copt544;
  Real copt546 = copt503 + copt62;
  Real copt547 = copt251 * copt546;
  Real copt549 = copt403 + copt502;
  Real copt550 = copt298 * copt549;
  Real copt551 = copt430 * copt515;
  Real copt552 = copt547 + copt550 + copt551;
  Real copt553 = copt552 * copt72;
  Real copt554 = copt528 + copt529 + copt530 + copt531 + copt532 + copt533 +
                 copt534 + copt535 + copt536 + copt545 + copt553;
  Real copt555 = copt439 * copt554;
  Real copt556 = ArcTan(copt527, copt555);
  Real copt559 = niexists(1);
  Real copt751 = -(copt49 * copt62);
  Real copt761 = -(copt77 * copt86);
  Real copt763 = -(copt251 * copt298);
  Real copt628 = copt508 * copt72;
  Real copt629 = copt611 + copt612 + copt613 + copt628;
  Real copt633 = copt632 * copt72;
  Real copt635 = copt47 + copt630;
  Real copt636 = copt635 * copt77;
  Real copt638 = copt404 * copt637;
  Real copt639 = copt633 + copt636 + copt638;
  Real copt640 = copt629 * copt639;
  Real copt641 = copt45 + copt631;
  Real copt642 = copt251 * copt641;
  Real copt644 = copt241 * copt643;
  Real copt646 = copt50 * copt645;
  Real copt647 = copt642 + copt644 + copt646;
  Real copt648 = copt492 * copt647;
  Real copt649 = -(copt251 * copt86);
  Real copt650 = copt241 * copt482;
  Real copt651 = copt436 * copt72;
  Real copt652 = copt298 * copt77;
  Real copt653 = copt649 + copt650 + copt651 + copt652;
  Real copt654 = -copt637;
  Real copt655 = copt654 + copt77;
  Real copt656 = copt241 * copt655;
  Real copt657 = copt637 + copt73;
  Real copt658 = copt251 * copt657;
  Real copt659 = copt407 * copt645;
  Real copt668 = copt656 + copt658 + copt659;
  Real copt669 = copt653 * copt668;
  Real copt670 = copt640 + copt648 + copt669;
  Real copt672 = copt251 * copt630 * copt86;
  Real copt673 = -(copt298 * copt630 * copt77);
  Real copt674 = copt251 * copt45 * copt637;
  Real copt675 = -(copt251 * copt62 * copt637);
  Real copt676 = -(copt298 * copt45 * copt637);
  Real copt677 = copt298 * copt49 * copt637;
  Real copt678 = copt632 * copt86;
  Real copt681 = copt417 + copt630;
  Real copt682 = copt681 * copt77;
  Real copt683 = copt508 * copt637;
  Real copt688 = copt678 + copt682 + copt683;
  Real copt690 = copt241 * copt688;
  Real copt691 = copt544 * copt645;
  Real copt693 = copt62 + copt631;
  Real copt694 = copt251 * copt693;
  Real copt700 = copt298 * copt643;
  Real copt701 = copt430 * copt645;
  Real copt702 = copt694 + copt700 + copt701;
  Real copt705 = copt702 * copt72;
  Real copt706 = copt528 + copt529 + copt672 + copt673 + copt674 + copt675 +
                 copt676 + copt677 + copt690 + copt691 + copt705;
  Real copt709 = copt415 * copt706;
  Real copt710 = ArcTan(copt670, copt709);
  Real copt713 = niexists(0);
  Real copt572 = copt571 * copt72;
  Real copt574 = copt573 * copt86;
  Real copt576 = copt418 * copt575;
  Real copt577 = copt572 + copt574 + copt576;
  Real copt578 = copt484 * copt577;
  Real copt579 = copt241 * copt571;
  Real copt580 = copt298 * copt573;
  Real copt582 = copt418 * copt581;
  Real copt583 = copt579 + copt580 + copt582;
  Real copt584 = copt492 * copt583;
  Real copt585 = -copt575;
  Real copt586 = copt585 + copt86;
  Real copt587 = copt241 * copt586;
  Real copt588 = copt575 + copt73;
  Real copt589 = copt298 * copt588;
  Real copt590 = copt421 * copt581;
  Real copt591 = copt587 + copt589 + copt590;
  Real copt592 = copt498 * copt591;
  Real copt593 = copt578 + copt584 + copt592;
  Real copt594 = copt251 * copt45 * copt86;
  Real copt595 = -(copt298 * copt45 * copt77);
  Real copt596 = -(copt251 * copt568 * copt86);
  Real copt597 = copt298 * copt568 * copt77;
  Real copt598 = -(copt251 * copt45 * copt575);
  Real copt599 = copt251 * copt575 * copt62;
  Real copt600 = copt298 * copt45 * copt575;
  Real copt601 = -(copt298 * copt49 * copt575);
  Real copt605 = copt571 * copt77;
  Real copt606 = copt403 + copt568;
  Real copt607 = copt606 * copt86;
  Real copt608 = copt430 * copt575;
  Real copt609 = copt605 + copt607 + copt608;
  Real copt610 = copt241 * copt609;
  Real copt614 = copt611 + copt612 + copt613;
  Real copt615 = copt581 * copt614;
  Real copt616 = copt49 + copt569;
  Real copt617 = copt298 * copt616;
  Real copt618 = copt417 + copt568;
  Real copt619 = copt251 * copt618;
  Real copt620 = copt508 * copt581;
  Real copt621 = copt617 + copt619 + copt620;
  Real copt622 = copt621 * copt72;
  Real copt623 = copt594 + copt595 + copt596 + copt597 + copt598 + copt599 +
                 copt600 + copt601 + copt610 + copt615 + copt622;
  Real copt624 = -(copt428 * copt623);
  Real copt625 = ArcTan(copt593, copt624);
  Real copt626 = niexists(2);
  Real copt801 = copt414 * copt441;
  Real copt802 = copt427 * copt467;
  Real copt806 = copt466 + copt801 + copt802;
  Real copt807 = 1 / copt806;
  Real copt808 = 1 / copt501;
  Real copt810 = copt427 * copt57;
  Real copt791 = copt42 * copt464;
  Real copt811 = copt791 + copt810;
  Real copt812 = Power(copt811, 2);
  Real copt815 = copt414 * copt42;
  Real copt788 = copt464 * copt57;
  Real copt816 = copt788 + copt815;
  Real copt824 = Power(copt816, 2);
  Real copt742 = -(copt72 * copt77);
  Real copt755 = copt45 * copt508;
  Real copt760 = copt72 * copt86;
  Real copt762 = copt241 * copt298;
  Real copt764 = copt447 + copt449 + copt451 + copt454 + copt742 + copt751 +
                 copt755 + copt760 + copt761 + copt762 + copt763;
  Real copt826 = copt42 * copt764;
  Real copt767 = copt241 * copt251;
  Real copt768 = copt430 * copt45;
  Real copt769 = copt432 * copt72;
  Real copt770 = copt462 + copt470 + copt472 + copt474 + copt751 + copt761 +
                 copt763 + copt767 + copt768 + copt769;
  Real copt827 = -(copt57 * copt770);
  Real copt828 = copt826 + copt827;
  Real copt831 = Power(copt828, 2);
  out(0)       = copt334;
  out(1)       = copt335 * copt398 *
           (copt207 * copt366 + copt311 * copt387 + copt357 * copt68);
  out(2) = copt397;
  out(3) = -(copt416 * copt429 * copt440 * copt478 * copt501 *
             (copt556 * copt558 * copt559 * copt560 * copt561 +
              copt567 * (copt441 * copt560 * copt625 * copt626 +
                         copt467 * copt561 * copt710 * copt713)));
  out(4) = copt334 * copt416 * copt429 * copt440 * copt730 * copt741 *
           (copt556 * copt557 * copt559 * copt560 * copt561 *
                (-(copt42 * copt764) + copt57 * copt770) +
            copt567 * (-(copt561 * copt57 * copt710 * copt713 *
                         (copt42 * copt452 + copt788)) +
                       copt42 * copt560 * copt625 * copt626 *
                           (copt475 * copt57 + copt791)));
  out(5) = copt807 * copt808 *
           (-(copt429 * copt625 * copt626 * copt812) -
            copt416 * copt710 * copt713 * copt824 -
            copt440 * copt556 * copt559 * copt831);

  return out;
}

#endif  // hylc_strain_II